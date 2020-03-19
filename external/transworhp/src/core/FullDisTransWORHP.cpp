#include "FullDisTransWORHP.h"

#include "conversion.h"
#include "xmlio.h"

#include "butcher.h"
#include "TWspline.h"
#include "TWfolder.h"
#include "TWproblem.h"
#include "twstatus.h"

#include "../base/vectortools.h"
#include "../base/defines.h"
#include "../base/exception.h"


#include <iomanip>
#include <sstream>
#include <algorithm>
#include <numeric>

#ifdef WIN32
#include <windows.h>
#endif

#undef max
#undef min

namespace tw {

using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::stringstream;

FullDisTransWorhp::FullDisTransWorhp(TransWorhpProblem *ph, TWparameter *twparam)
	: TransWorhp(ph,twparam),
	tempmatrixinitialised(false)
	{

	transworhp_type = TransWORHP_type::fullDiscretization;
}


FullDisTransWorhp::~FullDisTransWorhp() {

	if (tempmatrixinitialised) {
		FreeWorhpMatrix(&DFtemp1);
		FreeWorhpMatrix(&DGtemp1);
		FreeWorhpMatrix(&DFtemp2);
		FreeWorhpMatrix(&DGtemp2);
	}

}




void FullDisTransWorhp::GetSpace(int delta1, int delta2) {

	Delta1 = delta1;
	Delta2 = delta2;

	DS_obj.useStructure(phase->obj_structure(DS_obj));
	DS_ode.useStructure(phase->ode_structure(DS_ode));
	DS_rand.useStructure(phase->rand_structure(DS_rand));
	DS_neben.useStructure(phase->neben_structure(DS_neben));

	//if (n_integral)
	{
		DS_integral.useStructure(phase->integral_structure(DS_integral));
		if (DS_integral.useStructure()) {

			for (int k=0; k<n_integral; k++) {
				for (int i=0; i<twdiscretization->stuetzstellen(n_dis); i++) {
					for (int j=0; j<n_ode + n_ctrl; j++) {
						if (DS_integral.check(k,j)) {
							DS_obj(0, i * (n_ode+n_ctrl) + j);
						}
					}
				}

				for (int i=0; i<n_param; i++) {
					if (DS_integral.check(k,n_ode + n_ctrl + i)) {
						DS_obj(0, n_dis * (n_ode+n_ctrl) + i);
					}
				}
			}
		}
		DS_integral.finish();
	}

	DS_obj.finish();
	DS_ode.finish();
	DS_rand.finish();
	DS_neben.finish();

	printStructure();
}


double FullDisTransWorhp::x(int dis, int ode) const {

	if (dis<0) return 0;
	return X[(n_ode + n_ctrl) * twdiscretization->punkte() * dis + ode];
}
double &FullDisTransWorhp::setx(int dis, int ode)  {

	return X[(n_ode + n_ctrl) * twdiscretization->punkte() * dis + ode];
}
double FullDisTransWorhp::u(int dis, int ctrl) const {

	return X[(n_ode + n_ctrl) * twdiscretization->punkte() * dis + n_ode + ctrl];
}


double FullDisTransWorhp::p(int param) const {

	return X[(n_ode + n_ctrl) * twdiscretization->stuetzstellen(n_dis) + param];
}


double FullDisTransWorhp::x__(int dis, int ode) const {

	return X[(n_ode + n_ctrl) * dis + ode];
}

double FullDisTransWorhp::u__(int dis, int ctrl) const {

	return X[(n_ode + n_ctrl) * dis + n_ode + ctrl];
}




int FullDisTransWorhp::x_index(int dis, int ode) const {
	return (n_ode + n_ctrl) * twdiscretization->punkte() * dis + ode;
}

int FullDisTransWorhp::x_index__(int dis, int ode) const {
	return (n_ode + n_ctrl) * dis + ode;
}

int FullDisTransWorhp::u_index(int dis, int ctrl) const {
	return (n_ode + n_ctrl) * twdiscretization->punkte() * dis + n_ode + ctrl;
}

int FullDisTransWorhp::u_index__(int dis, int ctrl) const {
	return (n_ode + n_ctrl) * dis + n_ode + ctrl;
}

int FullDisTransWorhp::p_index(int param) const {
	return (n_ode + n_ctrl) * twdiscretization->stuetzstellen(n_dis) + param;
}

int FullDisTransWorhp::x_indexode(int ode) const {
	return ode;
}

int FullDisTransWorhp::u_indexode(int ctrl) const {
	return n_ode + ctrl;
}

int FullDisTransWorhp::p_indexode(int param) const {
	return n_ode + n_ctrl + param;
}


void FullDisTransWorhp::Connect(const OptVar &o, const Params &p) {

	//cout << "TW::CONN" << endl;

	Infty = p.Infty;

	X = &o.X[Delta1];
	X_low = &o.XL[Delta1];
	X_upp = &o.XU[Delta1];

	G = &o.G[Delta2 + 0];
	G_low = &o.GL[Delta2 + 0];
	G_upp = &o.GU[Delta2 + 0];

	Lambda = &o.Lambda[Delta1];
	Mu = &o.Mu[Delta2];

	for (int i=0; i<n_var; i++) {
		Lambda[i] = 0.0;
	}
	for (int i=0; i<n_con; i++) {
		Mu[i] = 0.0;
	}

	std::fill(tmp_ode_1.begin(),tmp_ode_1.end(),0.0);
	std::fill(tmp_ctrl_1.begin(),tmp_ctrl_1.end(),0.0);

	ZEN = &o.P[0];
	phase->zen_init(ZEN);

	for (int dis=0; dis<twdiscretization->stuetzstellen(n_dis); dis++) {

		phase->x_init(tmp_ode_1.data(),dis,twdiscretization->stuetzstellen(n_dis));

		const int index1 = (n_ode + n_ctrl) * dis;
		for (int ode = 0; ode < n_ode; ode++) {
			X[index1+ode] = tmp_ode_1[ode];
		}

		phase->u_init(tmp_ctrl_1.data(),dis,twdiscretization->stuetzstellen(n_dis));

		const int index2 = (n_ode + n_ctrl) * dis + n_ode;
		for (int ctrl = 0; ctrl < n_ctrl; ctrl++) {
			X[index2+ctrl] = tmp_ctrl_1[ctrl];
		}
	}

	phase->p_init(&X[twdiscretization->stuetzstellen(n_dis)*(n_ctrl+n_ode)]);

	phase->init();
}


void FullDisTransWorhp::Boundary() {

	std::fill(tmp_ode_1.begin(),tmp_ode_1.end(),-Infty);
	std::fill(tmp_ode_2.begin(),tmp_ode_2.end(),+Infty);

	std::fill(tmp_ctrl_1.begin(),tmp_ctrl_1.end(),-Infty);
	std::fill(tmp_ctrl_2.begin(),tmp_ctrl_2.end(),+Infty);

	phase->x_boundary(tmp_ode_1.data(),tmp_ode_2.data());
	phase->u_boundary(tmp_ctrl_1.data(),tmp_ctrl_2.data());

	for (int dis=0; dis<twdiscretization->stuetzstellen(n_dis); dis++) {

		for (int ode=0; ode<n_ode; ode++) {

			const int index = (n_ode + n_ctrl) * dis + ode;
			X_low[index] = tmp_ode_1[ode];
			X_upp[index] = tmp_ode_2[ode];
		}

		for (int ctrl=0; ctrl<n_ctrl; ctrl++) {

			const int index = (n_ode + n_ctrl) * dis + n_ode + ctrl;
			X_low[index] = tmp_ctrl_1[ctrl];
			X_upp[index] = tmp_ctrl_2[ctrl];
		}
	}

	for (int i = 0; i < n_param; i++) {
		X_low[twdiscretization->stuetzstellen(n_dis)*(n_ctrl+n_ode)+i] = -Infty;
		X_upp[twdiscretization->stuetzstellen(n_dis)*(n_ctrl+n_ode)+i] = +Infty;
	}
	phase->p_boundary( &X_low[twdiscretization->stuetzstellen(n_dis)*(n_ctrl+n_ode)] , &X_upp[twdiscretization->stuetzstellen(n_dis)*(n_ctrl+n_ode)] );

	phase->var_boundary(X_low, X_upp);

	// damit die Box-Schranken durch den Startwert nicht verletzt sind
	boxConToInitGuess();

	for (int i = 0; i < twdiscretization->stufen * (n_dis-1)*n_ode; i++) {
		G_low[i] = 0.0;
		G_upp[i] = 0.0;
	}

	if (n_rand) {
		std::fill(tmp_rand_1.begin(),tmp_rand_1.end(),0.0);
		std::fill(tmp_rand_2.begin(),tmp_rand_2.end(),0.0);

		phase->rand_boundary(tmp_rand_1.data(), tmp_rand_2.data());

		for (int i = 0; i < n_rand; i++) {
			const int index = twdiscretization->stufen * (n_dis-1) * n_ode + i;
			G_low[index] = tmp_rand_1[i];
			G_upp[index] = tmp_rand_2[i];
		}
	}

	if (n_neben) {
		std::fill(tmp_neben_1.begin(), tmp_neben_1.end(), -Infty);
		std::fill(tmp_neben_2.begin(), tmp_neben_2.end(), +Infty);

		phase->neben_boundary(tmp_neben_1.data(), tmp_neben_2.data());

		for (int i = 0; i < twdiscretization->stuetzstellen(n_dis); i++) {
			for (int j = 0; j < n_neben; j++) {
				const int index = twdiscretization->stufen * (n_dis-1) * n_ode + n_rand + i*n_neben + j;
				G_low[index] = tmp_neben_1[j];
				G_upp[index] = tmp_neben_2[j];
			}
		}
	}
}


void FullDisTransWorhp::RechteSeite(double *G1, double t, int dis) {

	double *X1 = &X[(n_ode+n_ctrl)* dis];
	double *U1 = &X[(n_ode+n_ctrl)* dis + n_ode];

	double *P = &X[(n_ode+n_ctrl)* twdiscretization->stuetzstellen(n_dis)];

	double DX1[1000];

	twcount_calls.ode.f += 1;

	phase->ode(DX1,t, //T[dis],
	    X1,  // Zustand X
	    U1, // Steuerung U
	    P);  // Parameter P



	for (int j=0; j<n_ode; j++) {
		G1[j] =  /*(T[dis+1] - T[dis]) *  */ (DX1[j]);
	}
}


void FullDisTransWorhp::Euler(double* G1, int dis) {

	double *X1 = &X[(n_ode+n_ctrl)* twdiscretization->punkte()*dis];
	double *U1 = &X[(n_ode+n_ctrl)* twdiscretization->punkte()*dis + n_ode];

	//double *X2 = &X[(n_ode+n_ctrl)* twdiscretization->punkte()*(dis+1)];
	//double *U2 = &X[(n_ode+n_ctrl)* twdiscretization->punkte()*(dis+1) + n_ode];

	double *P = &X[(n_ode+n_ctrl)* twdiscretization->stuetzstellen(n_dis)];

	double DX1[1000];//, DX2[1000];

	twcount_calls.ode.f +=2;
	phase->ode(DX1,T[dis],
	    X1,  // Zustand X
	    U1, // Steuerung U
	    P);  // Parameter P
	//ode(DX2,T[dis+1],
	 //   X2,  // Zustand X
	 //   U2, // Steuerung U
	 //   P); // Parameter P


	for (int j=0; j<n_ode; j++) {
		G1[j] = x(dis+1,j) - x(dis,j) - (T[dis+1] - T[dis]) * DX1[j];
	}
}


void FullDisTransWorhp::Trapez(double* G1, int dis) {

	double *X1 = &X[(n_ode+n_ctrl)* twdiscretization->punkte()*dis];
	double *U1 = &X[(n_ode+n_ctrl)* twdiscretization->punkte()*dis + n_ode];

	double *X2 = &X[(n_ode+n_ctrl)* twdiscretization->punkte()*(dis+1)];
	double *U2 = &X[(n_ode+n_ctrl)* twdiscretization->punkte()*(dis+1) + n_ode];

	double *P = &X[(n_ode+n_ctrl)* twdiscretization->stuetzstellen(n_dis)];

	double DX1[1000], DX2[1000];

	twcount_calls.ode.f +=2;
	phase->ode(DX1,T[dis],
	    X1,  // Zustand X
	    U1, // Steuerung U
	    P);  // Parameter P
	phase->ode(DX2,T[dis+1],
	    X2,  // Zustand X
	    U2, // Steuerung U
	    P); // Parameter P


	const double h = (T[dis+1] - T[dis]);
	for (int j=0; j<n_ode; j++) {
		G1[j] = x(dis+1,j) - x(dis,j) -  h * (DX1[j]+DX2[j])/2;
	}
}


void FullDisTransWorhp::HermiteSimpson(double* G1, double *G2, int dis) {

	double *X1 = &X[(n_ode+n_ctrl)* twdiscretization->punkte()*dis];
	double *U1 = &X[(n_ode+n_ctrl)* twdiscretization->punkte()*dis + n_ode];

	double *X2 = &X[(n_ode+n_ctrl)* twdiscretization->punkte()*(dis+1)];
	double *U2 = &X[(n_ode+n_ctrl)* twdiscretization->punkte()*(dis+1) + n_ode];

	double *X12 = &X[(n_ode+n_ctrl)* (twdiscretization->punkte()*dis+1)];
	double *U12 = &X[(n_ode+n_ctrl)* (twdiscretization->punkte()*dis+1) + n_ode];

	double *P = &X[(n_ode+n_ctrl)* twdiscretization->stuetzstellen(n_dis)];

	double DX1[1000], DX2[1000],DX12[1000];

	twcount_calls.ode.f +=3;
	phase->ode(DX1,T[dis],
	    X1,  // Zustand X
	    U1, // Steuerung U
	    P);  // Parameter P
	phase->ode(DX2,T[dis+1],
	    X2,  // Zustand X
	    U2, // Steuerung U
	    P);  // Parameter P
	phase->ode(DX12,(T[dis]+T[dis+1])*.5,
	    X12,  // Zustand X
	    U12, // Steuerung U
	    P);  // Parameter P



	// Hermite
	for (int j=0; j<n_ode; j++) {
		G1[j] = X[(n_ode + n_ctrl) * (2 * dis+1) + j]
		        - .5*(x(dis+1,j) + x(dis,j))
		        - .125*(T[dis+1] - T[dis]) * (DX1[j]-DX2[j]);

	}

	// Simpson
	for (int j=0; j<n_ode; j++) {
		G2[j] = x(dis+1,j) - x(dis,j) - (T[dis+1] - T[dis])/6. * (DX1[j]+DX2[j]+4.*DX12[j]);


	}
}

void FullDisTransWorhp::Lobatto(double* G1, double *G2, double *G3, int dis) {

	double *X1 = &X[(n_ode+n_ctrl)* twdiscretization->punkte()*dis];
	double *U1 = &X[(n_ode+n_ctrl)* twdiscretization->punkte()*dis + n_ode];

	double *X2 = &X[(n_ode+n_ctrl)* twdiscretization->punkte()*(dis+1)];
	double *U2 = &X[(n_ode+n_ctrl)* twdiscretization->punkte()*(dis+1) + n_ode];

	double *XM = &X[(n_ode+n_ctrl)* (twdiscretization->punkte()*dis+1)];
	double *UM = &X[(n_ode+n_ctrl)* (twdiscretization->punkte()*dis+1) + n_ode];

	double *XP = &X[(n_ode+n_ctrl)* (twdiscretization->punkte()*dis+2)];
	double *UP = &X[(n_ode+n_ctrl)* (twdiscretization->punkte()*dis+2) + n_ode];

	double *P = &X[(n_ode+n_ctrl)* twdiscretization->stuetzstellen(n_dis)];

	double DX1[1000], DX2[1000],DXM[1000], DXP[1000];

	twcount_calls.ode.f +=3;
	phase->ode(DX1, T[dis],
	    X1,  // Zustand X
	    U1, // Steuerung U
	    P);  // Parameter P
	phase->ode(DX2, T[dis+1],
	    X2,  // Zustand X
	    U2, // Steuerung U
	    P);  // Parameter P
	phase->ode(DXM, T[dis] + (T[dis+1]-T[dis])*(5.-sqrt(5.))/10.,
	    XM,  // Zustand X
	    UM, // Steuerung U
	    P);  // Parameter P
	phase->ode(DXP, T[dis] + (T[dis+1]-T[dis])*(5.+sqrt(5.))/10.,
	    XP,  // Zustand X
	    UP, // Steuerung U
	    P);  // Parameter P


	for (int j=0; j<n_ode; j++) {

		double x1 = x(dis,j);
		double x2 = X[(n_ode + n_ctrl) * (3 * dis+1) + j];
		double x3 = X[(n_ode + n_ctrl) * (3 * dis+2) + j];
		double x4 = x(dis+1,j);
		double h = T[dis+1] - T[dis];
		double f1 = DX1[j];
		double f2 = DXM[j];
		double f3 = DXP[j];
		double f4 = DX2[j];

		G1[j] =	(-1./6.+(5./6.)*sqrt(5.))*x1 - (31./21.)*x2 + (23./14.-(5./6.)*sqrt(5.))*x3
				+ h * ( (-1./20.+(43./420.)*sqrt(5.))*f1 + ((2./105.)*sqrt(5.)-1./30.)*f4);

		G2[j] = (13./6.+(5./6.)*sqrt(5.))*x1 + (-23./14.-(5./6.)*sqrt(5.))*x2 - (11./21.)*x3
				+ h * ( (7./30.+(3./35.)*sqrt(5.))*f1 + (1./60.+(1./420.)*sqrt(5.))*f4);

		G3[j] = x4 - x1 - h/12. * (f1 + 5.*f2 + 5.*f3 + f4);
	}



}

void FullDisTransWorhp::Constraints2(double *GG, int /*DGflag*/) {

	if (twdiscretization->type==TWdiscretizationType::Trapez) {

		for (int dis=0; dis<n_dis-1; dis++) {
			double *G1 = &GG[dis*n_ode];
			Trapez(G1, dis);
		}
	} else if (twdiscretization->type==TWdiscretizationType::HermiteSimpson) {

		for (int dis=0; dis<n_dis-1; dis++) {
			double *G1 = &GG[2*dis*n_ode];
			double *G2 = &GG[(2*dis+1)*n_ode];

			HermiteSimpson(G1,G2, dis);
		}

	} else if (twdiscretization->type==TWdiscretizationType::Lobatto) {

		for (int dis=0; dis<n_dis-1; dis++) {
			double *G1 = &GG[3*dis*n_ode];
			double *G2 = &GG[(3*dis+1)*n_ode];
			double *G3 = &GG[(3*dis+2)*n_ode];

			Lobatto(G1,G2,G3, dis);

		}
	} else if (twdiscretization->type==TWdiscretizationType::Euler) {

		for (int dis=0; dis<n_dis-1; dis++) {
			double *G1 = &GG[dis*n_ode];
			Euler(G1, dis);
		}
	}

	if (n_rand) {
		twcount_calls.rand.f++;
		phase->rand( &GG[twdiscretization->stufen * (n_dis-1) * n_ode] );
	}

	if (n_neben) {
		for (int dis=0; dis<twdiscretization->stuetzstellen(n_dis); dis++) {
			double *G1 = &GG[twdiscretization->stufen * (n_dis-1) * n_ode + n_rand + dis*n_neben];

			double *X1 = &X[(n_ode+n_ctrl) * dis];
			double *U1 = &X[(n_ode+n_ctrl) * dis + n_ode];
			double *P = &X[(n_ode+n_ctrl)* twdiscretization->stuetzstellen(n_dis)];
			twcount_calls.neben.f++;

			if (twdiscretization->type==TWdiscretizationType::Trapez || twdiscretization->type==TWdiscretizationType::Euler) {
				phase->neben(G1,T[dis],X1,U1,P);
			}
			else if (twdiscretization->type==TWdiscretizationType::HermiteSimpson) {

				double tt;

				if (dis%2==0) {
					tt = T[dis/2];
				} else {
					tt = (T[(dis-1)/2] + T[(dis-1)/2+1])/2.;
				}

				phase->neben(G1,tt,X1,U1,P);
			}


		}
	}

}

double FullDisTransWorhp::Objective(double ScaleObj) {

	twcount_calls.obj.f++;

	Lagrange();
	double ret = phase->obj();
	for (int i=0; i<n_integral; i++)
		ret += lagrange_weight[i] * lagrange_integral[n_integral*(twdiscretization->stuetzstellen(n_dis)-1)+i];
	//cout << ret << endl;
	return ret * ScaleObj;
}



int FullDisTransWorhp::DG_structure(const TWfolder *f, WorhpMatrix *DG, int offset) {

	DG_start = offset;
	int ind=offset;
	double val;
	double mult;

	// Ableitungen nach x und u fuer alle diskreten Punkte
	for (int k=0; k<n_dis; k++) {

		// und Hilfspunkte
		for (int kk=0; kk<twdiscretization->punkte(); kk++) {

			if (kk>0 && k==n_dis-1) break;

			// Ableitungen nach x
			for (int l=0; l<n_ode; l++) {

				int colindex = (twdiscretization->punkte()*k+kk)*(n_ode + n_ctrl) +  l + 1;

				// Ableitungen des ODE-Systems
				if (k>0 && kk==0) // Block A
					for (int i=0; i<twdiscretization->stufen; i++) {
						for (int j=0; j<n_ode; j++) {
							if (DS_ode.odecheck(j,l, twdiscretization->punkte(), i, val,mult,twdiscretization->type)) {
								if (DG) DG->col[ind] = colindex + Delta1;
								if (DG) DG->row[ind] = (k-1)*n_ode*twdiscretization->stufen + n_ode*i + j +1 + Delta2;
								if (DG) DG->val[ind] = 0;
								ind++;
							}
						}
					}

				if (k<n_dis-1) // Block B
					for (int i=0; i<twdiscretization->stufen; i++) {
						for (int j=0; j<n_ode; j++) {
							if (DS_ode.odecheck(j,l,kk,i,val,mult,twdiscretization->type)) {
								if (DG) DG->col[ind] = colindex + Delta1;
								if (DG) DG->row[ind] = k*n_ode*twdiscretization->stufen + n_ode*i + j + 1 + Delta2;
								if (DG) DG->val[ind] = 0;
								ind++;
							}
						}
					}

				//if (kk==0)
				{
					// Ableitungen nach Randwerten
					for (int j=0; j<n_rand; j++) {

						//int userindex = k*(n_ode + n_ctrl) +  l;

						if (DS_rand.check(j,colindex-1))
						{
							if (DG) DG->col[ind] = colindex + Delta1;
							if (DG) DG->row[ind] = ( twdiscretization->stuetzstellen(n_dis) - 1 )*n_ode + j + 1 + Delta2;
							if (DG) DG->val[ind] = 0;
							ind++;
						}
					}
				}
				if (k<n_dis-1 || kk==0) {
					// Ableitungen nach Nebenbedingungen
					for (int j=0; j<n_neben; j++) {
						if (DS_neben.check(j,l)) {
							if (DG) DG->col[ind] = colindex + Delta1;
							if (DG) DG->row[ind] = ( twdiscretization->stuetzstellen(n_dis) - 1 )*n_ode + n_rand + (twdiscretization->punkte()*k+kk)*n_neben + j + 1 + Delta2;
							//if (DG) DG->row[ind] = (n_dis-1)*n_ode + n_rand + k*n_neben + j + 1;
							if (DG) DG->val[ind] = 0;
							ind++;
						}
					}
				}

// Zentrale ABLEITUNGEN
if (f) // Struktur der Beschraenkungen im TWFolder
				f->DG_structure(ind,colindex + Delta1, DG);


			}


			// Ableitungen nach u
			for (int l=0; l<n_ctrl; l++) {

				int colindex = (twdiscretization->punkte()*k+kk)*(n_ode + n_ctrl) + n_ode +  l + 1;

				// Ableitungen des ODE-Systems
				if (k>0 && kk==0) // Block C
					for (int i=0; i<twdiscretization->stufen; i++) {
						for (int j=0; j<n_ode; j++) {
							if (DS_ode.odecheck(j,n_ode+l, twdiscretization->punkte(), i, val,mult,twdiscretization->type)) {
								//if (DS_ode.check(j,n_ode+l)) {
								if (DG) DG->col[ind] = colindex + Delta1;
								if (DG) DG->row[ind] = (k-1)*n_ode*twdiscretization->stufen+ n_ode*i + j + 1 + Delta2;
//((k-1)*n_ode + j )*stufen+i+ 1;
								if (DG) DG->val[ind] = 0;
								ind++;
							}
						}
					}

				if (k<n_dis-1) // Block D
					for (int i=0; i<twdiscretization->stufen; i++) {
						for (int j=0; j<n_ode; j++) {
							if (DS_ode.odecheck(j,n_ode+l, kk, i, val,mult,twdiscretization->type)) {
								//if (DS_ode.check(j,n_ode+l)) {
								if (DG) DG->col[ind] = colindex + Delta1;
								if (DG) DG->row[ind] = k*n_ode*twdiscretization->stufen+ n_ode*i + j + 1 + Delta2;
//(k*n_ode + j )*stufen+i+ 1;
								if (DG) DG->val[ind] = 0;
								ind++;
							}
						}
					}

				//if (kk==0)
				{
					// Ableitungen nach Randwerten
					for (int j=0; j<n_rand; j++) {

						// int colindex = (Disc.punkte()*k+kk)*(n_ode + n_ctrl) + n_ode +  l + 1;
						//int userindex = k*(n_ode + n_ctrl) + n_ode +  l;


						if (DS_rand.check(j,colindex-1))
						{
							if (DG) DG->col[ind] = colindex + Delta1;
							if (DG) DG->row[ind] = ( twdiscretization->stuetzstellen(n_dis) - 1 )*n_ode + j + 1 + Delta2;
							// (n_dis-1)*n_ode + j + 1;
							if (DG) DG->val[ind] = 0;
							ind++;
						}
					}
				}

				// Ableitungen nach Nebenbedingungen
				if (k<n_dis-1 || kk==0) {
					for (int j=0; j<n_neben; j++) {
						if (DS_neben.check(j,n_ode+l)) {
							if (DG) DG->col[ind] = colindex + Delta1;
							if (DG) DG->row[ind] = ( twdiscretization->stuetzstellen(n_dis) - 1 )*n_ode + n_rand + (twdiscretization->punkte()*k+kk)*n_neben + j + 1 + Delta2;
							//   if (DG) DG->row[ind] = (n_dis-1)*n_ode + n_rand + k*n_neben + j + 1;
							if (DG) DG->val[ind] = 0;
							ind++;
						}
					}
				}


// Zentrale ABLEITUNGEN
				if (f) f->DG_structure(ind,colindex + Delta1, DG);

			}

		}
	}


	// Ableitungen nach Parametern
	for (int l=0; l<n_param; l++) {

		int colindex = twdiscretization->stuetzstellen(n_dis)*(n_ode + n_ctrl) +  l + 1;

		for (int k=0; k<n_dis-1; k++) {

			for (int i=0; i<twdiscretization->stufen; i++) {
				for (int j=0; j<n_ode; j++) {

					if (DS_ode.check(j,n_ode+n_ctrl + l)) {
						if (DG) DG->col[ind] = colindex + Delta1;
						if (DG) DG->row[ind] = k*n_ode*twdiscretization->stufen + i*n_ode + j + 1+ Delta2; // ???????????
// k*n_ode*stufen+ n_ode*i + j + 1;
						if (DG) DG->val[ind] = 0;
						ind++;
					}
				}
			}
		}

		for (int j=0; j<n_rand; j++) {
			if (DS_rand.check(j,colindex-1)) {
				if (DG) DG->col[ind] = colindex + Delta1;
				if (DG) DG->row[ind] =  ( twdiscretization->stuetzstellen(n_dis) - 1 )*n_ode + j + 1 + Delta2;
				//(n_dis-1)*n_ode + j + 1;
				if (DG) DG->val[ind] = 0;
				ind++;
			}
		}

		for (int k=0; k<n_dis; k++) {
			for (int kk=0; kk<twdiscretization->punkte(); kk++) {
				if (k<n_dis-1 || kk==0) {

					for (int j=0; j<n_neben; j++) {
						if (DS_neben.check(j,n_ode+n_ctrl + l)) {
							if (DG) DG->col[ind] = colindex + Delta1;
							if (DG) DG->row[ind] = ( twdiscretization->stuetzstellen(n_dis) - 1 )*n_ode + n_rand + (twdiscretization->punkte()*k+kk)*n_neben + j + 1 + Delta2;
							//(n_dis-1)*n_ode + n_rand + j + 1;
							if (DG) DG->val[ind] = 0;
							ind++;
						}
					}
				}
			}
		}

// Zentrale ABLEITUNGEN
				if (f) f->DG_structure(ind,colindex + Delta1, DG);



	}


	/*if (DG) {
	for (int i=0;i<ind;i++)
	 print(*DG,i);
	}
	cout << "DG_structure: " << ind << endl;
	*/

	DG_nnz = ind - DG_start;
	return ind;


}



// Ableitung df_k/d(y_k,u_k,param) erzeugen
void FullDisTransWorhp::DG_diff_ode(double t, int dis, int active_index) {

	// in DS_ode ablegen
	// dis von 0 bis Disc.stuetzstellen(n_dis)

	DS_ode.activeindex = active_index;

	// Analytisch?
	double *X1 = &X[(n_ode+n_ctrl)* dis];
	double *U1 = &X[(n_ode+n_ctrl)* dis + n_ode];
	double *P = &X[(n_ode+n_ctrl)* twdiscretization->stuetzstellen(n_dis)];

	int ret  = phase->ode_diff(DS_ode, t, X1, U1, P);
	if (ret) {
		twcount_calls.ode.df++;
		return;
	}

	// oder numerisch
	for (int l=0; l<n_ode+n_ctrl; l++) {

		int var_index = dis*(n_ode + n_ctrl) +  l;

		// zentraler Differenzenquotient
		double v = X[var_index];
		X[var_index] += eps;
		RechteSeite(tmp_ode_1.data(),t,dis);
		X[var_index] = v - eps;
		RechteSeite(tmp_ode_2.data(),t,dis);
		X[var_index] = v;


	//	ScaleODE(tmp_ode_1);
	//	ScaleODE(tmp_ode_2);

		for (int j=0; j<n_ode; j++) {

			double *de = DS_ode.find(j,l);
			if (de) de[active_index] =  (tmp_ode_1[j]-tmp_ode_2[j]) / (2*eps);

			//			if (DS_ode.check(j,l))
			//	DS_ode.set(j,l,  (tmp_ode_1[j]-tmp_ode_2[j]) / (2*eps) );
			//dg[j][l] = (tmp1[j]-tmp2[j]) / (2*eps);
		}
	}
}




void FullDisTransWorhp::DG_diff_ode_p(double t, int l, int dis, int active_index) {

	// in DS_ode ablegen
	// dis von 0 bis Disc.stuetzstellen(n_dis)

	DS_ode.activeindex = active_index;

	double *X1 = &X[(n_ode+n_ctrl)* dis];
	double *U1 = &X[(n_ode+n_ctrl)* dis + n_ode];
	double *P = &X[(n_ode+n_ctrl)* twdiscretization->stuetzstellen(n_dis)];

	int ret = phase->ode_diff_p(DS_ode, t, X1, U1, P, n_ode+n_ctrl+l);
	if (ret) {
		twcount_calls.ode.dfp++;
		return;
	}

	int var_index = twdiscretization->stuetzstellen(n_dis)*(n_ode + n_ctrl) + l;

	double v = X[var_index];
	X[var_index] += eps;
	RechteSeite(tmp_ode_1.data(),t,dis);
	X[var_index] = v - eps;
	RechteSeite(tmp_ode_2.data(),t,dis);
	X[var_index] = v;

	//	ScaleODE(tmp_ode_1);
	//	ScaleODE(tmp_ode_2);

	for (int j=0; j<n_ode; j++) {

		double *de = DS_ode.find(j,n_ode + n_ctrl + l);
		if (de) de[active_index] = (tmp_ode_1[j]-tmp_ode_2[j]) / (2*eps);

		//if (DS_ode.check()) {
		//	DS_ode.set(j,n_ode + n_ctrl + l,  (tmp_ode_1[j]-tmp_ode_2[j]) / (2*eps) );
		//}
	}
}




void FullDisTransWorhp::DG_diff_rand(int var_index) {

	// dis von 0 bis Disc.stuetzstellen(n_dis)
	//cout << "DG_calculate " << dis<< endl;



	//int var_index = Disc.stuetzstellen(n_dis)*(n_ode + n_ctrl) + l;

	twcount_calls.rand.f+=2;

	double v = X[var_index];
	X[var_index] += eps;
	phase->rand(tmp_rand_1.data());
	X[var_index] = v - eps;
	phase->rand(tmp_rand_2.data());
	X[var_index] = v;

	for (int j=0; j<n_rand; j++) {

		double *de = DS_rand.find(j,var_index);
		if (de) de[0] =  (tmp_rand_1[j]-tmp_rand_2[j]) / (2*eps);

		// dr[j] = (tmp_rand_1[j]-tmp_rand_2[j]) / (2*eps);
	}
}



void FullDisTransWorhp::DG_diff_neben(double t, int dis) {

	// in DS_neben ablegen
	// dis von 0 bis Disc.stuetzstellen(n_dis)

	if (!n_neben) return;


	// Analytisch?
	double *X1 = &X[(n_ode+n_ctrl)* dis];
	double *U1 = &X[(n_ode+n_ctrl)* dis + n_ode];
	double *P = &X[(n_ode+n_ctrl)* twdiscretization->stuetzstellen(n_dis)];

	int ret = phase->neben_diff(DS_neben, t, X1, U1, P);
	if (ret) {
		twcount_calls.neben.df++;
		return;
	}

	// oder numerisch
	for (int l=0; l<n_ode+n_ctrl; l++) {

		int var_index = dis*(n_ode + n_ctrl) +  l ;

		// zentraler Differenzenquotient
		twcount_calls.neben.f+=2;
		double v = X[var_index];
		X[var_index] += eps;
		phase->neben(tmp_neben_1.data(),t,X1,U1,P);
		X[var_index] = v - eps;
		phase->neben(tmp_neben_2.data(),t,X1,U1,P);
		X[var_index] = v;
		for (int j=0; j<n_neben; j++) {

			double *de = DS_neben.find(j,l);
			if (de) de[0] =  (tmp_neben_1[j]-tmp_neben_2[j]) / (2*eps);


			//	if (s.check(j,l))
			//	s.set(j,l,(tmp_neben_1[j]-tmp_neben_2[j]) / (2*eps) );
		}
	}
}




void FullDisTransWorhp::DG_diff_neben_p(double t, int l, int dis) {

	// in DS_neben ablegen
	// dis von 0 bis Disc.stuetzstellen(n_dis)

	if (!n_neben) return;

	double *X1 = &X[(n_ode+n_ctrl)* dis];
	double *U1 = &X[(n_ode+n_ctrl)* dis + n_ode];
	double *P = &X[(n_ode+n_ctrl)* twdiscretization->stuetzstellen(n_dis)];

	int ret = phase->neben_diff_p(DS_neben, t, X1, U1, P, n_ode+n_ctrl+l);
	if (ret) {
		twcount_calls.neben.dfp++;
		return;
	}

	int var_index = twdiscretization->stuetzstellen(n_dis)*(n_ode + n_ctrl) + l;

	double v = X[var_index];
	X[var_index] += eps;
	phase->neben(tmp_neben_1.data(),t,X1,U1,P);
	X[var_index] = v - eps;
	phase->neben(tmp_neben_2.data(),t,X1,U1,P);
	X[var_index] = v;

	for (int j=0; j<n_neben; j++) {

		double *de = DS_neben.find(j,n_ode + n_ctrl + l);
		if (de) de[0] =  (tmp_neben_1[j]-tmp_neben_2[j]) / (2*eps);

//		dn[j] = (tmp_neben_1[j]-tmp_neben_2[j]) / (2*eps);
	}
}


void FullDisTransWorhp::DG_calculate(TWfolder *f, WorhpMatrix &DG) {

	int ind=0;

	double val;
	double mult;


	int use_auto_rand = phase->rand_diff(DS_rand);
	if (use_auto_rand) twcount_calls.rand.df++;


	for (int k=0; k<n_dis; k++) {

		if (twdiscretization->type==TWdiscretizationType::Trapez) {
			DG_diff_ode(T[k], twdiscretization->punkte()*k, 0);
			if (k<n_dis-1) DG_diff_ode(T[k+1],twdiscretization->punkte()*k+1, 1); //?

		} else if (twdiscretization->type==TWdiscretizationType::HermiteSimpson) {
			DG_diff_ode(T[k], twdiscretization->punkte()*k, 0);
			if (k<n_dis-1) {
				DG_diff_ode(.5 * (T[k] + T[k+1]),twdiscretization->punkte()*k+1, 1);
				DG_diff_ode(T[k+1],twdiscretization->punkte()*k+2, 2);
			}

		} else if (twdiscretization->type==TWdiscretizationType::Lobatto) {
			DG_diff_ode(T[k], twdiscretization->punkte()*k, 0);
			if (k<n_dis-1) {
				DG_diff_ode(T[k] + (T[k+1]-T[k]) * (5 - sqrt(5.))/10,twdiscretization->punkte()*k+1, 1);
				DG_diff_ode(T[k] + (T[k+1]-T[k]) * (5 + sqrt(5.))/10,twdiscretization->punkte()*k+2, 1);
				DG_diff_ode(T[k+1],twdiscretization->punkte()*k+3, 2);
			}

		} else if (twdiscretization->type==TWdiscretizationType::Euler) {
			DG_diff_ode(T[k], twdiscretization->punkte()*k, 0);
			//if (k<n_dis-1) DG_diff_ode(T[k+1],twdiscretization->punkte()*k+1, 1); //?
		}



		//if (k<5)
		/*{
			cout << "K=" << k<< endl;
			for (int z=0;z<3;z++) {
				for (int j=0;j<n_ode;j++) {
					for (int l=0;l<n_ode+n_ctrl+n_param;l++) {


						cout << setw(10) << dg_[z][j][l];
					}
					cout << endl;

				}
			}
		}*/


		for (int kk=0; kk<twdiscretization->punkte(); kk++) {

				if (kk>0 && k==n_dis-1) break;


			//virtual bool neben_diff(DiffStructure &s, double t, const double *x, const double *u, const double *p) {return false;}



			if (twdiscretization->type==TWdiscretizationType::Trapez) {

				if (kk==0) DG_diff_neben(T[k],twdiscretization->punkte()*k);
				// DG_calculate(dg_[0],T[k],  Disc.punkte()*k);
				//else if (kk==1 && k<n_dis-1)
				//   DG_diff_neben(DS_neben,T[k+1],Disc.punkte()*k+1);
				// DG_calculate(dg_[1],T[k+1],Disc.punkte()*k+1);

			} else if (twdiscretization->type==TWdiscretizationType::HermiteSimpson) {
				if (kk==0) DG_diff_neben(T[k],twdiscretization->punkte()*k);
				// DG_calculate(dg_[0],T[k],Disc.punkte()*k);
				else if (kk==1 && k<n_dis-1)
					DG_diff_neben(.5 * (T[k] + T[k+1]),twdiscretization->punkte()*k+1);
				//DG_calculate(dg_[1],.5 * (T[k] + T[k+1]),Disc.punkte()*k+1);
				//else
				//if (k<n_dis-1) DG_calculate(dg_[2],T[k+1],Disc.punkte()*k+2);

			} else if (twdiscretization->type==TWdiscretizationType::Lobatto) {
				if (kk==0) DG_diff_neben(T[k],twdiscretization->punkte()*k);
				// DG_calculate(dg_[0],T[k],Disc.punkte()*k);
				else if (kk==1 && k<n_dis-1)
					DG_diff_neben(T[k] + (T[k+1]-T[k])*(5.-sqrt(5.)/10),twdiscretization->punkte()*k+1);
				else if (kk==2 && k<n_dis-1)
					DG_diff_neben(T[k] + (T[k+1]-T[k])*(5.+sqrt(5.)/10),twdiscretization->punkte()*k+2);
				//DG_calculate(dg_[1],.5 * (T[k] + T[k+1]),Disc.punkte()*k+1);
				//else
				//if (k<n_dis-1) DG_calculate(dg_[2],T[k+1],Disc.punkte()*k+2);

			} else if (twdiscretization->type==TWdiscretizationType::Euler) {
				if (kk==0) DG_diff_neben(T[k],twdiscretization->punkte()*k);
			}



			for (int l=0; l<n_ode; l++) {

				int colindex = (twdiscretization->punkte()*k+kk)*(n_ode + n_ctrl) +  l + 1;

				if (k>0 && kk==0) // Block A
					for (int i=0; i<twdiscretization->stufen; i++) {
						for (int j=0; j<n_ode; j++) {
							if (DS_ode.odecheck(j,l, twdiscretization->punkte(), i, val,mult,twdiscretization->type)) {
								//if (DG) DG->col[ind] = colindex;
								//if (DG) DG->row[ind] = (k-1)*n_ode*Disc.stufen + n_ode*i + j +1;
//((k-1)*n_ode + j )*stufen+i +1;

								DS_ode.activeindex=0;
								DG.val[DG_start + ind] = val + (T[k]-T[k-1])*mult* DS_ode.get(j,l);
								//dg_[0/*Disc.punkte()*/][j][l]; // vom vorherigen!
								ind++;
							}
						}
					}

				if (k<n_dis-1) // Block B
					for (int i=0; i<twdiscretization->stufen; i++) {
						for (int j=0; j<n_ode; j++) {
							if (DS_ode.odecheck(j,l,kk,i,val,mult,twdiscretization->type)) {
								//if (DG) DG->col[ind] = colindex;
								//if (DG) DG->row[ind] = k*n_ode*Disc.stufen + n_ode*i + j + 1;
//(k*n_ode + j )*stufen+i + 1;
								DS_ode.activeindex=kk;
								DG.val[DG_start + ind] = val + (T[k+1]-T[k])*mult* DS_ode.get(j,l);
								//dg_[kk][j][l];
								ind++;
							}
						}
					}



				//if (kk==0)
				{


					if (n_rand && !use_auto_rand) DG_diff_rand(colindex-1);

					for (int j=0; j<n_rand; j++) {

						double *de = DS_rand.find(j,colindex-1);
						if (de) {
							//if (DS_rand.check(j,colindex-1)) {
							// if (DG) DG->col[ind] = colindex;
							// if (DG) DG->row[ind] = (n_dis-1)*n_ode + j + 1;
							// DG.val[ind] = DS_rand.get(j,colindex-1);
							DG.val[DG_start + ind] = de[0];
							ind++;
						}
					}
				}

				if (k<n_dis-1 || kk==0) {

					//int varindex = (Disc.punkte()*k+kk)*(n_ode + n_ctrl) + l;
					//double tt = T[k];
					//if (Disc.mode==1 && kk==1)
//						tt =.5 * (T[k] + T[k+1]);
					//				DG_diff_neben(dn,tt,  colindex-1,(Disc.punkte()*k+kk));

					// DS_neben


					for (int j=0; j<n_neben; j++) {
						if (DS_neben.check(j,l)) {
							// if (DG) DG->col[ind] = colindex;
							// if (DG) DG->row[ind] = (n_dis-1)*n_ode + n_rand + k*n_neben + j + 1;
							DG.val[DG_start + ind] = DS_neben.get(j,l);
							ind++;
						}
					}

				}


	if (f) f->DG_calculate(DG_start, ind, colindex + Delta1, DG);


			}



			for (int l=0; l<n_ctrl; l++) {

				int colindex = (twdiscretization->punkte()*k+kk)*(n_ode + n_ctrl) + n_ode +  l + 1;

				if (k>0 && kk==0) // Block C
					for (int i=0; i<twdiscretization->stufen; i++) {
						for (int j=0; j<n_ode; j++) {
							if (DS_ode.odecheck(j,n_ode+l, twdiscretization->punkte(), i, val,mult,twdiscretization->type)) {
								//if (DS_ode.check(j,n_ode+l)) {
								//  if (DG) DG->col[ind] = colindex;
								//  if (DG) DG->row[ind] = (k-1)*n_ode*Disc.stufen+ n_ode*i + j + 1;
//((k-1)*n_ode + j )*stufen+i+ 1;
								DS_ode.activeindex=0;
								DG.val[DG_start + ind] = val + (T[k]-T[k-1])*mult*DS_ode.get(j,n_ode+l);
								//dg_[0/*Disc.punkte()*/][j][n_ode+l];
								ind++;
							}
						}
					}

				if (k<n_dis-1) // Block D
					for (int i=0; i<twdiscretization->stufen; i++) {
						for (int j=0; j<n_ode; j++) {
							if (DS_ode.odecheck(j,n_ode+l, kk, i, val,mult,twdiscretization->type)) {
								//if (DS_ode.check(j,n_ode+l)) {
								//  if (DG) DG->col[ind] = colindex;
								//  if (DG) DG->row[ind] = k*n_ode*Disc.stufen+ n_ode*i + j + 1;
//(k*n_ode + j )*stufen+i+ 1;
								DS_ode.activeindex=kk;
								DG.val[DG_start + ind] = val + (T[k+1]-T[k])*mult*DS_ode.get(j,n_ode+l);
								//dg_[kk][j][n_ode+l];
								ind++;
							}
						}
					}


				//if (kk==0)
				{

					if (n_rand && !use_auto_rand) DG_diff_rand(colindex-1);

					for (int j=0; j<n_rand; j++) {

						double *de = DS_rand.find(j,colindex-1);
						if (de) {


							//int userindex = k*(n_ode + n_ctrl) + n_ode +  l;

							//if (DS_rand.check(j,colindex-1)) {
							// if (DG) DG->col[ind] = colindex;
							// if (DG) DG->row[ind] = (n_dis-1)*n_ode + j + 1;
							//	DG.val[ind] = DS_rand.get(j,colindex-1);
							DG.val[DG_start + ind] = de[0];
							//	DG.val[ind] = dr[j];
							ind++;
						}
					}
				}

				if (k<n_dis-1 || kk==0) {
					//int colindex = (twdiscretization->punkte()*k+kk)*(n_ode + n_ctrl) + n_ode +  l + 1;


					//double tt = T[k];
					//if (Disc.mode==1 && kk==1)
					//	tt =.5 * (T[k] + T[k+1]);
					//DG_diff_neben(dn_,tt,  colindex-1,(Disc.punkte()*k+kk));


					for (int j=0; j<n_neben; j++) {
						if (DS_neben.check(j,n_ode+l)) {
							// if (DG) DG->col[ind] = colindex;
							// if (DG) DG->row[ind] = (n_dis-1)*n_ode + n_rand + k*n_neben + j + 1;
							//				DG.val[ind] = dn_[j];
							DG.val[DG_start + ind] = DS_neben.get(j,n_ode+l);

							ind++;
						}
					}

				}


				if (f) f->DG_calculate(DG_start,ind, colindex + Delta1, DG);

			}



		}
	}

	for (int l=0; l<n_param; l++) {

		int colindex = twdiscretization->stuetzstellen(n_dis)*(n_ode + n_ctrl) +  l + 1;

		for (int k=0; k<n_dis-1; k++) {


			//double dgp_[3][100];
			//int varindex =  l;


			if (twdiscretization->type==TWdiscretizationType::Trapez) {
				DG_diff_ode_p(T[k],   l, k,   0);
				DG_diff_ode_p(T[k+1], l, k+1, 1); //?

			} else if (twdiscretization->type==TWdiscretizationType::HermiteSimpson) {
				DG_diff_ode_p(      T[k],           l, twdiscretization->punkte()*k,   0);
				DG_diff_ode_p(.5 * (T[k] + T[k+1]), l, twdiscretization->punkte()*k+1, 1);
				DG_diff_ode_p(             T[k+1] , l, twdiscretization->punkte()*k+2, 2);

			} else if (twdiscretization->type==TWdiscretizationType::Lobatto) {
				DG_diff_ode_p(      T[k],           l, twdiscretization->punkte()*k,   0);
				DG_diff_ode_p(T[k] + (T[k+1] - T[k])*(5-sqrt(5.))/10., l, twdiscretization->punkte()*k+1, 1);
				DG_diff_ode_p(T[k] + (T[k+1] - T[k])*(5+sqrt(5.))/10., l, twdiscretization->punkte()*k+2, 2);
				DG_diff_ode_p(             T[k+1] , l, twdiscretization->punkte()*k+3, 3);

			} else if (twdiscretization->type==TWdiscretizationType::Euler) {
				DG_diff_ode_p(T[k],   l, k,   0);
			}

//if (k<5)
			/*{
			 cout << "K=" << k<< " L="<< l << endl;
			 for (int z=0;z<3;z++) {
			  for (int j=0;j<n_ode;j++) {
			    cout << setw(10) << dgp_[z][j];

			   cout << endl;

			  }
			 }
			}*/



			for (int i=0; i<twdiscretization->stufen; i++) {
				for (int j=0; j<n_ode; j++) {

					if (DS_ode.odecheck(j,n_ode+n_ctrl+l, 0, i, val,mult,twdiscretization->type)) {
						//(j,n_ode+n_ctrl + l)) {
						// if (DG) DG->col[ind] = colindex;
						// if (DG) DG->row[ind] = k*n_ode*Disc.stufen + i*n_ode + j + 1; // ???????????
// k*n_ode*stufen+ n_ode*i + j + 1;

						double *de = DS_ode.find(j,n_ode+n_ctrl+l);
						if (!de) throw Exception("ODE");


						if (twdiscretization->type==TWdiscretizationType::Trapez) {
							DG.val[DG_start + ind] = (T[k+1]-T[k])*mult*(de[0]+de[1]);

						} else if (twdiscretization->type==TWdiscretizationType::HermiteSimpson) {
							if (i==0)
								DG.val[DG_start + ind] = (T[k+1]-T[k])*mult*(de[0]-de[2]);
							if (i==1)
								DG.val[DG_start + ind] = (T[k+1]-T[k])*mult*(de[0]+4.*de[1]+de[2]);

						} else if (twdiscretization->type==TWdiscretizationType::Lobatto) {
							cout << "TODO Lobatto" << endl;
						//	if (i==0)
						//		DG.val[DG_start + ind] = (T[k+1]-T[k])*mult*(de[0]-de[2]);
						//	if (i==1)
						//		DG.val[DG_start + ind] = (T[k+1]-T[k])*mult*(de[0]+4.*de[1]+de[2]);

						} else if (twdiscretization->type==TWdiscretizationType::Euler) {
							DG.val[DG_start + ind] = (T[k+1]-T[k])*mult*(de[0]);
						}



						ind++;
					}
				}
			}
		}


		if (n_rand && !use_auto_rand) DG_diff_rand(colindex-1);


		for (int j=0; j<n_rand; j++) {

			double *de = DS_rand.find(j,colindex-1);
			if (de) {


				//if (DS_rand.check(j,colindex-1)) {
				// if (DG) DG->col[ind] = colindex;
				// if (DG) DG->row[ind] = (n_dis-1)*n_ode + j + 1;
				//DG.val[ind] = dr[j];
				//				DG.val[ind] = DS_rand.get(j,colindex-1);
				DG.val[DG_start + ind] = de[0];
				ind++;
			}
		}


		for (int k=0; k<n_dis; k++) {
			for (int kk=0; kk<twdiscretization->punkte(); kk++) {
				if (k<n_dis-1 || kk==0) {


					//int varindex = (Disc.punkte()*k+kk)*(n_ode + n_ctrl) + l;

					double tt = T[k];
					if (twdiscretization->type==TWdiscretizationType::HermiteSimpson && kk==1)
						tt =.5 * (T[k] + T[k+1]);

					if (twdiscretization->type==TWdiscretizationType::Lobatto && kk==1)
						tt = T[k] + (T[k+1] - T[k]) * (5.-sqrt(5.))/10.;

					if (twdiscretization->type==TWdiscretizationType::Lobatto && kk==2)
						tt = T[k] + (T[k+1] - T[k]) * (5.+sqrt(5.))/10.;

					//..... colindex!!! falsch!
					//DG_diff_neben_p(tt,  colindex-1,(twdiscretization->punkte()*k+kk));
DG_diff_neben_p(tt,  l,(twdiscretization->punkte()*k+kk));

					for (int j=0; j<n_neben; j++) {

						double *de = DS_neben.find(j,n_ode+n_ctrl + l);

						if (de) {
							DG.val[DG_start + ind] = de[0];
							ind++;
						}
						//if (DS_neben.check(j%n_neben,n_ode+n_ctrl + l)) {
						// if (DG) DG->col[ind] = colindex;
						// if (DG) DG->row[ind] = (n_dis-1)*n_ode + n_rand + j + 1;
						//DG.val[ind] = dn_[j];
						//ind++;
						//}

					}
				}
			}
		}
		//}

	if (f)	f->DG_calculate(DG_start, ind, colindex + Delta1, DG);


	}


}


int FullDisTransWorhp::HM_structure_ohne_Diag(int hessianstructure, WorhpMatrix &/*DF*/, WorhpMatrix &/*DG*/, WorhpMatrix *HM, int offset) {

	HM_start = offset;

	//cout << "HM_structure_ohne_Diag" << HM_start << endl;
	int ind=offset;


	if (hessianstructure==0) { // Diagonal


	}

	else if (hessianstructure==1) { // Full

		for (int y = 0; y < n_var; ++y) {
			for (int x = y+1; x < n_var; ++x) {

				int use = 1;

				if (use) {
					if (HM) {
						HM->row[ind] = x + 1 + Delta1;
						HM->col[ind] = y + 1 + Delta1;
					}
					ind++;
				}
			}
		}
	}

	else if (hessianstructure==2) { // Blockode2, erweitert auf benachbarte Bloecke


		for (int y = 0; y < n_var; ++y) {
			for (int x = y+1; x < n_var; ++x) {

				int use = 0;

				int xblock = x / ((n_ode+n_ctrl));
				int yblock = y / ((n_ode+n_ctrl));

				if (std::abs(xblock-yblock)<2) use=1; // Benachbarte Bloecke verbinden

				if (x >= (n_ode+n_ctrl)*twdiscretization->stuetzstellen(n_dis) ) use=1; // mit freien Par. verbinden

				if (use) {
					if (HM) {
						HM->row[ind] = x + 1 + Delta1;
						HM->col[ind] = y + 1 + Delta1;
					}
					ind++;
				}
			}
		}
	}

	else if (hessianstructure==3) { // Blockode

		int blocksize = (n_ode+n_ctrl);

		for (int y = 0; y < n_var; ++y) {
			for (int x = y+1; x < n_var; ++x) {

				int use = 0;

				int xblock = x / blocksize;
				int yblock = y / blocksize;

				if (std::abs(xblock-yblock)<1) use=1; // Id. Bloecke verbinden

				if (x >= (n_ode+n_ctrl)*twdiscretization->stuetzstellen(n_dis) ) use=1; // mit freien Par. verbinden

				if (use) {
					if (HM) {
						HM->row[ind] = x + 1 + Delta1;
						HM->col[ind] = y + 1 + Delta1;
					}
					ind++;
				}
			}
		}
	}

	HM_nnz_ohne_diag = ind - HM_start;

	return ind;
}


void FullDisTransWorhp::tempmatrixinit() {

	if (!tempmatrixinitialised) {
		tempmatrixinitialised=true;

		int Delta1_ = Delta1;
		int Delta2_ = Delta2;
		Delta1 = 0;
		Delta2 = 0;


		ZeroWorhpMatrix(&DFtemp1);
		DFtemp1.nnz = DF_structure();
		DFtemp1.nRow = n_var;
		DFtemp1.nCol = 1;
		DFtemp1.kind = 5;
		InitWorhpMatrix(&DFtemp1, (char*)"DFtemp1", 0,0,0);
		DF_structure(&DFtemp1);

		ZeroWorhpMatrix(&DFtemp2);
		DFtemp2.nnz = DFtemp1.nnz;
		DFtemp2.nRow = n_var;
		DFtemp2.nCol = 1;
		DFtemp2.kind = 5;
		InitWorhpMatrix(&DFtemp2, (char*)"DFtemp2", 0,0,0);
		DF_structure(&DFtemp2);

		ZeroWorhpMatrix(&DGtemp1);
		DGtemp1.nnz = DG_structure(nullptr); // nur fuer Hessematrix, d.h. 2. Ableitung = 0
		DGtemp1.nRow = n_con;
		DGtemp1.nCol = n_var;
		DGtemp1.kind = 1;
		InitWorhpMatrix(&DGtemp1, (char*)"DGtemp1", 0,0,0);
		DG_structure(nullptr,&DGtemp1);

		ZeroWorhpMatrix(&DGtemp2);
		DGtemp2.nnz = DGtemp1.nnz;
		DGtemp2.nRow = n_con;
		DGtemp2.nCol = n_var;
		DGtemp2.kind = 1;
		InitWorhpMatrix(&DGtemp2, (char*)"DGtemp2", 0,0,0);
		DG_structure(nullptr,&DGtemp2);


		Delta1 = Delta1_;
		Delta2 = Delta2_;
	}
}


void FullDisTransWorhp::HM_calculate5(WorhpMatrix &/*DF*/, WorhpMatrix &/*DG*/, WorhpMatrix &/*HM*/, double /*ScaleObj*/, double */*Mu*/) {
}

void FullDisTransWorhp::HM_calculate4(WorhpMatrix &/*DF*/, WorhpMatrix &/*DG*/, WorhpMatrix &/*HM*/, double /*ScaleObj*/, double */*Mu*/) {

	// "(values=2) Calculating 1st derivative of DF + mu*DG .. faster",
exit(1);
/*
	WorhpMatrix DFtemp1,DGtemp1,DFtemp2,DGtemp2;

//	double mudg[1000];

	ZeroWorhpMatrix0(&DFtemp1);
	DFtemp1.nnz = DF.nnz;
	DFtemp1.nRow = worhp_o->n;
	DFtemp1.nCol = 1;
	DFtemp1.kind = 5;
	InitWorhpMatrix0(&DFtemp1, (char*)"DFtemp1", 0,0,0);
	DF_structure(&DFtemp1);

	ZeroWorhpMatrix0(&DFtemp2);
	DFtemp2.nnz = DF.nnz;
	DFtemp2.nRow = worhp_o->n;
	DFtemp2.nCol = 1;
	DFtemp2.kind = 5;
	InitWorhpMatrix0(&DFtemp2, (char*)"DFtemp2", 0,0,0);
	DF_structure(&DFtemp2);

	ZeroWorhpMatrix0(&DGtemp1);
	DGtemp1.nnz = DG.nnz;
	DGtemp1.nRow = worhp_o->m;
	DGtemp1.nCol = worhp_o->n;
	DGtemp1.kind = 1;
	InitWorhpMatrix0(&DGtemp1, (char*)"DGtemp1", 0,0,0);
	DG_structure(&DGtemp1);

	ZeroWorhpMatrix0(&DGtemp2);
	DGtemp2.nnz = DG.nnz;
	DGtemp2.nRow = worhp_o->m;
	DGtemp2.nCol = worhp_o->n;
	DGtemp2.kind = 1;
	InitWorhpMatrix0(&DGtemp2, (char*)"DGtemp2", 0,0,0);
	DG_structure(&DGtemp2);

	double eps = 1e-6;

	int ind = 0;

	for (int x=1; x<=worhp_o->n; x++) {


		double tmp_x = worhp_o->X[x-1];

		worhp_o->X[x-1] = tmp_x + eps;

		DF_calculate(DFtemp1,ScaleObj);
		DG_calculate(DGtemp1);

		worhp_o->X[x-1] = tmp_x - eps;

		DF_calculate(DFtemp2,ScaleObj);
		DG_calculate(DGtemp2);

		worhp_o->X[x-1] = tmp_x;


		while (HM.col[ind] == x) {

			int y = HM.row[ind];




			double lll = 0; //

			for (int i=0; i<DFtemp1.nnz; i++) {

				if (DFtemp1.row[i] == y) {
					lll = (DFtemp1.val[i] - DFtemp2.val[i]) / (2*eps);
					break;
				}
			}

			for (int i=0; i<DGtemp1.nnz; i++) {


				if (DGtemp1.col[i] == y) {
					double zz = (DGtemp1.val[i] - DGtemp2.val[i])/(2*eps) * Mu[DGtemp1.row[i]-1];
					lll += zz;
					//			cout << "(" << setw(2) << x << "/" << setw(2) << y << ")   " << setw(10) << i
					//			<< setw(20) << DGtemp1.val[i] << setw(20) << DGtemp2.val[i] << setw(20) << zz << endl;

				}

			}

			HM.val[ind] = lll;

			ind++;
		}



		{
			int y = x;




			double lll = 0; //

			for (int i=0; i<DFtemp1.nnz; i++) {

				if (DFtemp1.row[i] == y) {
					lll = (DFtemp1.val[i] - DFtemp2.val[i]) / (2*eps);
					break;
				}
			}

			for (int i=0; i<DGtemp1.nnz; i++) {


				if (DGtemp1.col[i] == y) {
					double zz = (DGtemp1.val[i] - DGtemp2.val[i])/(2*eps) * Mu[DGtemp1.row[i]-1];
					lll += zz;
					//			cout << "(" << setw(2) << x << "/" << setw(2) << y << ")   " << setw(10) << i
					//			<< setw(20) << DGtemp1.val[i] << setw(20) << DGtemp2.val[i] << setw(20) << zz << endl;

				}

			}



			HM.val[HM.nnz - worhp_o->n + x-1] = lll;


		}



	}
*/
}

void FullDisTransWorhp::HM_calculate3(WorhpMatrix &/*DF*/, WorhpMatrix &/*DG*/, WorhpMatrix &HM, double ScaleObj, double *Mu) {

	//cout << "HM_calc3" << endl;
	for (int i=0;i<n_var;i++) {
		cout << std::setw(20) << Lambda[i];
		if ((i%5)==4) cout << endl;
	}
	for (int i=0;i<n_con;i++) {
		cout << std::setw(20) << Mu[i];
		if ((i%5)==4) cout << endl;
	}

	// "(values=2) Calculating 1st derivative of DF + mu*DG .. faster",

	tempmatrixinit();


	 int DF_start_ = DF_start;
	 int DG_start_ = DG_start;

	 DF_start = DG_start = 0;

	double eps = 1e-6;

	int ind = 0;

	for (int x=1; x<=n_var; x++) {


		double tmp_x = X[x-1];

		X[x-1] = tmp_x + eps;

		DF_calculate(DFtemp1,ScaleObj);
		DG_calculate(nullptr,DGtemp1);

		X[x-1] = tmp_x - eps;

		DF_calculate(DFtemp2,ScaleObj);
		DG_calculate(nullptr,DGtemp2);

		X[x-1] = tmp_x;


		while (HM.col[HM_start + ind] == x + Delta1) {

			int y = HM.row[HM_start + ind] - Delta1;

			double lll = 0; //

			for (int i=0; i<DFtemp1.nnz; i++) {

				if (DFtemp1.row[i] == y) {
					lll = (DFtemp1.val[i] - DFtemp2.val[i]) / (2*eps);

					///cout << "HM_DF " << i << " " << lll << endl;
					break;
				}
			}

			for (int i=0; i<DGtemp1.nnz; i++) {


				if (DGtemp1.col[i] == y) {
					double zz = (DGtemp1.val[i] - DGtemp2.val[i])/(2*eps) * Mu[DGtemp1.row[i] -1];
					lll += zz;

					///cout << "HM_DG " << i << " " << zz << endl;
					//			cout << "(" << setw(2) << x << "/" << setw(2) << y << ")   " << setw(10) << i
					//			<< setw(20) << DGtemp1.val[i] << setw(20) << DGtemp2.val[i] << setw(20) << zz << endl;

				}

			}

			HM.val[HM_start + ind] = lll;

			ind++;
		}


		 {
			int y = x;




			double lll = 0; //

			for (int i=0; i<DFtemp1.nnz; i++) {

				if (DFtemp1.row[i] == y) {
					lll = (DFtemp1.val[i] - DFtemp2.val[i]) / (2*eps);

					///	cout << "HM_DF: " << i << " " << lll << endl;

					break;
				}
			}

			for (int i=0; i<DGtemp1.nnz; i++) {


				if (DGtemp1.col[i] == y) {
					double zz = (DGtemp1.val[i] - DGtemp2.val[i])/(2*eps) * Mu[DGtemp1.row[i]-1];
					lll += zz;






		//	if ( Mu[DGtemp1.row[i]-1] )
		//	cout << "Mu3 "<< Mu[DGtemp1.row[i]-1] << setw(20) << Mu[DGtemp1.row[i]-1] << endl;


						/**cout << "HM_DG: " << i << " " << zz << "   " << Mu[DGtemp1.row[i]+Delta2-1]
						<< "::    " << DGtemp1.val[i] << " - " << DGtemp2.val[i] << " = "
						<< DGtemp1.val[i] - DGtemp2.val[i] << endl;
				*/
					//			cout << "(" << setw(2) << x << "/" << setw(2) << y << ")   " << setw(10) << i
					//			<< setw(20) << DGtemp1.val[i] << setw(20) << DGtemp2.val[i] << setw(20) << zz << endl;

				}

			}




			HM.val[/*HM.nnz - n_var + x-1*/   +  HM_start_diag + x-1  ] = lll;


		}



	}


	MatrixToMATLAB(HM,"test.m");

	DF_start = DF_start_;
	DG_start = DG_start_;




}

void FullDisTransWorhp::HM_calculate2(WorhpMatrix &/*DF*/, WorhpMatrix &/*DG*/, WorhpMatrix &HM, double ScaleObj, double *Mu) {

	// "(values=1) Calculating 1st derivative of DF + mu*DG",

	tempmatrixinit();

	double eps = 1e-6;

	for (int ind=0; ind<HM_nnz_ohne_diag; ind++) {

		int x = HM.col[HM_start + ind];
		int y = HM.row[HM_start + ind];

		double tmp_y = X[y-1];

		X[y-1] = tmp_y + eps;

		DF_calculate(DFtemp1,ScaleObj);
		DG_calculate(nullptr,DGtemp1);

		X[y-1] = tmp_y - eps;

		DF_calculate(DFtemp2,ScaleObj);
		DG_calculate(nullptr,DGtemp2);

		X[y-1] = tmp_y;

		double lll = 0; //

		for (int i=0; i<DFtemp1.nnz; i++) {

			if (DFtemp1.row[i] == x) {
				lll = (DFtemp1.val[i] - DFtemp2.val[i]) / (2*eps);
				break;
			}
		}
		//cout << "lll: " << x << " " << y << " " << DFtemp1.val[x-1] <<" " << DFtemp1.val[x-1]<< " " << lll<< endl;

		for (int i=0; i<DGtemp1.nnz; i++) {


			if (DGtemp1.col[i] == x) {
				double zz = (DGtemp1.val[i] - DGtemp2.val[i])/(2*eps) * Mu[DGtemp1.row[i]-1];
				lll += zz;
				//			cout << "(" << setw(2) << x << "/" << setw(2) << y << ")   " << setw(10) << i
				//			<< setw(20) << DGtemp1.val[i] << setw(20) << DGtemp2.val[i] << setw(20) << zz << endl;

			}

		}

		//cout << "(" << setw(2) << x << "/" << setw(2) << y << ")   ";
		//cout << "------ " << tmp_y << " " << " ------> " << lll << endl;

		HM.val[HM_start_diag + ind] = lll;




	}

}



void FullDisTransWorhp::HM_calculate1(WorhpMatrix &/*DF*/, WorhpMatrix &/*DG*/, WorhpMatrix &HM, double ScaleObj, double *Mu) {
	//cout << "HM_calc1" << endl;
	// "(values=0) Calculating 2nd derivative of F + mu*G",

	for (int ind=0; ind<HM_nnz_ohne_diag; ind++) {

		double eps = 1e-6;
		int x = HM.col[HM_start + ind] - Delta1;
		int y = HM.row[HM_start + ind] - Delta1;

		double ddf;

	/*	if (x==y) {
			double tmp_x = X[x-1];

			X[x-1] = tmp_x + eps;
			double f1 = Objective(ScaleObj);
			Constraints2(tmp_gg_1);

			X[x-1] = tmp_x - eps;
			double f3 = Objective(ScaleObj);
			Constraints2(tmp_gg_3);

			X[x-1] = tmp_x;
			double f2 = Objective(ScaleObj);
			Constraints2(tmp_gg_2);

			ddf = (f1-2*f2+f3) / (eps*eps);

			for (int i=0; i<n_con; i++) {
				tmp_gg_5[i] = (tmp_gg_1[i]-2*tmp_gg_2[i]+tmp_gg_3[i]) / (eps*eps);
			}

		}

		else {*/
			double tmp_x = X[x-1];
			double tmp_y = X[y-1];

			X[x-1] = tmp_x + eps;
			X[y-1] = tmp_y + eps;
			double f1 = Objective(ScaleObj);
			Constraints2(tmp_gg_1.data());

			X[x-1] = tmp_x + eps;
			X[y-1] = tmp_y - eps;
			double f2 = Objective(ScaleObj);
			Constraints2(tmp_gg_2.data());

			X[x-1] = tmp_x - eps;
			X[y-1] = tmp_y + eps;
			double f3 = Objective(ScaleObj);
			Constraints2(tmp_gg_3.data());

			X[x-1] = tmp_x - eps;
			X[y-1] = tmp_y - eps;
			double f4 = Objective(ScaleObj);
			Constraints2(tmp_gg_4.data());

			X[x-1] = tmp_x;
			X[y-1] = tmp_y;

			ddf = (f1-f2-f3+f4) / (4*eps*eps);

			for (int i=0; i<n_con; i++) {
				tmp_gg_5[i] = (tmp_gg_1[i]-tmp_gg_2[i]-tmp_gg_3[i]+tmp_gg_4[i]) / (4*eps*eps);
			}
		//}


		double lxx = ddf;

		for (int i=0; i<n_con; i++) {
			lxx += Mu[i] * tmp_gg_5[i];
		}


		HM.val[HM_start + ind] = lxx;

	}





	for (int ind=0; ind<n_var; ind++) {

		double eps = 1e-6;
		int x = HM.col[HM_start_diag + ind] - Delta1;
		//int y = HM.row[HM_start_diag + ind] - Delta1;

		double ddf;

		//if (x==y) {
			double tmp_x = X[x-1];

			X[x-1] = tmp_x + eps;
			double f1 = Objective(ScaleObj);
			Constraints2(tmp_gg_1.data());

			X[x-1] = tmp_x - eps;
			double f3 = Objective(ScaleObj);
			Constraints2(tmp_gg_3.data());

			X[x-1] = tmp_x;
			double f2 = Objective(ScaleObj);
			Constraints2(tmp_gg_2.data());

			ddf = (f1-2*f2+f3) / (eps*eps);

			for (int i=0; i<n_con; i++) {
				tmp_gg_5[i] = (tmp_gg_1[i]-2*tmp_gg_2[i]+tmp_gg_3[i]) / (eps*eps);
			}

		/*}

		else {
			double tmp_x = X[x-1];
			double tmp_y = X[y-1];

			X[x-1] = tmp_x + eps;
			X[y-1] = tmp_y + eps;
			double f1 = Objective(ScaleObj);
			Constraints2(tmp_gg_1);

			X[x-1] = tmp_x + eps;
			X[y-1] = tmp_y - eps;
			double f2 = Objective(ScaleObj);
			Constraints2(tmp_gg_2);

			X[x-1] = tmp_x - eps;
			X[y-1] = tmp_y + eps;
			double f3 = Objective(ScaleObj);
			Constraints2(tmp_gg_3);

			X[x-1] = tmp_x - eps;
			X[y-1] = tmp_y - eps;
			double f4 = Objective(ScaleObj);
			Constraints2(tmp_gg_4);

			X[x-1] = tmp_x;
			X[y-1] = tmp_y;

			ddf = (f1-f2-f3+f4) / (4*eps*eps);

			for (int i=0; i<n_con; i++) {
				tmp_gg_5[i] = (tmp_gg_1[i]-tmp_gg_2[i]-tmp_gg_3[i]+tmp_gg_4[i]) / (4*eps*eps);
			}
		}*/


		double lxx = ddf;

	//	cout << setw(20) << lxx << endl;
		for (int i=0; i<n_con; i++) {
			lxx += Mu[i] * tmp_gg_5[i];

		//	if (Mu[i])
		//	cout << "Mu"<< i << setw(20) << Mu[i] << setw(20) << tmp_gg_5[i] << endl;
		}


		HM.val[HM_start_diag + ind] = lxx;

	}



}


void FullDisTransWorhp::Lagrange() {

	if (!n_integral) return;

	// fuer Trapezregel:
	std::fill_n(lagrange_integral.begin(),n_integral,0.0);

	double *P = &X[(n_ode+n_ctrl)* twdiscretization->stuetzstellen(n_dis)];

	if (twdiscretization->type==TWdiscretizationType::Euler) { // wie bei Trapez!?
		for (int dis=0; dis<n_dis-1; dis++) {

			double *X1 = &X[(n_ode+n_ctrl)* twdiscretization->punkte()*dis];
			double *U1 = &X[(n_ode+n_ctrl)* twdiscretization->punkte()*dis + n_ode];

			double *X2 = &X[(n_ode+n_ctrl)* twdiscretization->punkte()*(dis+1)];
			double *U2 = &X[(n_ode+n_ctrl)* twdiscretization->punkte()*(dis+1) + n_ode];

			phase->integral(tmp_integral_1.data(), T[dis], X1, U1, P);
			phase->integral(tmp_integral_2.data(), T[dis+1], X2, U2, P);

			for (int k=0; k<n_integral; k++) {
				lagrange_integral[n_integral*(dis+1) + k ] = lagrange_integral[n_integral*(dis) + k ]
				        + (tmp_integral_1[k]+tmp_integral_2[k])/2 * (T[dis+1]-T[dis]);
			}
		}

	}

	else if (twdiscretization->type==TWdiscretizationType::Trapez) {
		for (int dis=0; dis<n_dis-1; dis++) {

			double *X1 = &X[(n_ode+n_ctrl)* twdiscretization->punkte()*dis];
			double *U1 = &X[(n_ode+n_ctrl)* twdiscretization->punkte()*dis + n_ode];

			double *X2 = &X[(n_ode+n_ctrl)* twdiscretization->punkte()*(dis+1)];
			double *U2 = &X[(n_ode+n_ctrl)* twdiscretization->punkte()*(dis+1) + n_ode];

			phase->integral(tmp_integral_1.data(), T[dis], X1, U1, P);
			phase->integral(tmp_integral_2.data(), T[dis+1], X2, U2, P);

			for (int k=0; k<n_integral; k++) {
				lagrange_integral[n_integral*(dis+1) + k ] = lagrange_integral[n_integral*(dis) + k ]
				        + (tmp_integral_1[k]+tmp_integral_2[k])/2 * (T[dis+1]-T[dis]);
			}
		}

	}

	else if (twdiscretization->type==TWdiscretizationType::HermiteSimpson) {
		for (int dis=0; dis<n_dis-1; dis++) { // nur Intervalle abgehen!

			double *X1 = &X[(n_ode+n_ctrl)* twdiscretization->punkte()*dis];
			double *U1 = &X[(n_ode+n_ctrl)* twdiscretization->punkte()*dis + n_ode];

			double *X2 = &X[(n_ode+n_ctrl)* twdiscretization->punkte()*(dis+1)];
			double *U2 = &X[(n_ode+n_ctrl)* twdiscretization->punkte()*(dis+1) + n_ode];

			double *X12 = &X[(n_ode+n_ctrl)* (twdiscretization->punkte()*dis+1)];
			double *U12 = &X[(n_ode+n_ctrl)* (twdiscretization->punkte()*dis+1) + n_ode];

			double t1 = T[dis];
			double t2 = T[dis+1];
			double t12 = (t1+t2)/2;

			phase->integral(tmp_integral_1.data(), t1, X1, U1, P);
			phase->integral(tmp_integral_12.data(), t12, X12, U12, P);
			phase->integral(tmp_integral_2.data(), t2, X2, U2, P);

			for (int k=0; k<n_integral; k++) {
				lagrange_integral[n_integral*(2*dis+2) + k] = lagrange_integral[n_integral*(2*dis)+k]
				        + (tmp_integral_1[k] + 4*tmp_integral_12[k] + tmp_integral_2[k])/6 * (t2-t1);

				lagrange_integral[n_integral*(2*dis+1) + k] = lagrange_integral[n_integral*(2*dis)+k]
				        + (5*tmp_integral_1[k] + 8*tmp_integral_12[k] - tmp_integral_2[k])/24 * (t2-t1);

			}
		}

	}

	else if (twdiscretization->type==TWdiscretizationType::Lobatto) {

		cout << "TODO 2 Lobatto" << endl;
	/*	for (int dis=0; dis<n_dis-1; dis++) { // nur Intervalle abgehen!

			double *X1 = &X[(n_ode+n_ctrl)* twdiscretization->punkte()*dis];
			double *U1 = &X[(n_ode+n_ctrl)* twdiscretization->punkte()*dis + n_ode];

			double *X2 = &X[(n_ode+n_ctrl)* twdiscretization->punkte()*(dis+1)];
			double *U2 = &X[(n_ode+n_ctrl)* twdiscretization->punkte()*(dis+1) + n_ode];

			double *X12 = &X[(n_ode+n_ctrl)* (twdiscretization->punkte()*dis+1)];
			double *U12 = &X[(n_ode+n_ctrl)* (twdiscretization->punkte()*dis+1) + n_ode];

			double t1 = T[dis];
			double t2 = T[dis+1];
			double t12 = (t1+t2)/2;

			integral(tmp_integral_1.data(), t1, X1, U1, P);
			integral(tmp_integral_12.data(), t12, X12, U12, P);
			integral(tmp_integral_2.data(), t2, X2, U2, P);

			for (int k=0; k<n_integral; k++) {
				lagrange_integral[n_integral*(2*dis+2) + k] = lagrange_integral[n_integral*(2*dis)+k]
				        + (tmp_integral_1[k] + 4*tmp_integral_12[k] + tmp_integral_2[k])/6 * (t2-t1);

				lagrange_integral[n_integral*(2*dis+1) + k] = lagrange_integral[n_integral*(2*dis)+k]
				        + (5*tmp_integral_1[k] + 8*tmp_integral_12[k] - tmp_integral_2[k])/24 * (t2-t1);

			}
		}*/

	}

}



// evtl. bei Integrationsverfahren anpassen...
int FullDisTransWorhp::DF_structure(WorhpMatrix *DF, int offset) {

	//cout << "DF_structure" << offset << endl;

	DF_start = offset;
	int ind = offset;
	for (int k = 0; k < n_var; k++) {

		if (DS_obj.check(0,k)) {

			if (DF) {
				DF->row[ind] = k+1 + Delta1;
				DF->val[ind] = 0;
			}
			ind++;
		}
	}


	//for (int i=0;i<ind;i++)
	// if (DF) print(*DF,i);

	//cout << "DF_size = " <<  ind << endl;
	DF_nnz = ind - DF_start;
	return ind;
}


void FullDisTransWorhp::DF_calculate(WorhpMatrix &DF, double ScaleObj) {


	if (DS_integral.useStructure()) {
		for (int k=0; k<n_integral; k++) {
			for (int j=0; j<n_ode + n_ctrl; j++) {
				if (DS_integral.check(k,j)) {
					for (int i=0; i<twdiscretization->stuetzstellen(n_dis); i++) {
						DS_obj(0, i * (n_ode+n_ctrl) + j) = 0;
					}
				}
			}
		}
	}

	/*
		for (int i=0;i<n_var;i++) {
			if (DS_obj.check(0,i)) DS_obj(0, i) = 0;
		}*/

	// Ableitung von obj zur Verfuegung?
	bool ret = phase->obj_diff(DS_obj);


	if (ret && n_integral > 0) {

		int i=0;
		double *X1 = &X[(n_ode+n_ctrl)* i];
		double *U1 = &X[(n_ode+n_ctrl)* i + n_ode];

		double *P = &X[(n_ode+n_ctrl)* twdiscretization->stuetzstellen(n_dis)];

		// auch Ableitung von integral zur Verfuegung?
		bool ret2 = phase->integral_diff(DS_integral,T[i],X1,U1,P);
		if (ret2) {

			// fuer Trapez!!!

			///lagrange[k] += (f1[k]+f2[k])/2 * (T[i]-T[i-1]);

			if (twdiscretization->type==TWdiscretizationType::Euler) { // genauso hochintegrieren wie bei Trapez!?
				for (int i=0; i<n_dis; i++) {

					double *X1 = &X[(n_ode+n_ctrl)* i];
					double *U1 = &X[(n_ode+n_ctrl)* i + n_ode];

					double *P = &X[(n_ode+n_ctrl)* twdiscretization->stuetzstellen(n_dis)];

					phase->integral_diff(DS_integral,T[i],X1,U1,P);

					for (int k=0; k<n_integral; k++) {
						for (int j=0; j<n_ode + n_ctrl; j++) {
							if (DS_integral.check(k,j)) {

								if (i>0)
									DS_obj(0, i * (n_ode+n_ctrl) + j) += (T[i]-T[i-1])/2 * DS_integral(k,j)*lagrange_weight[k];
								if (i<n_dis-1)
									DS_obj(0, i * (n_ode+n_ctrl) + j) += (T[i+1]-T[i])/2 * DS_integral(k,j)*lagrange_weight[k];
							}
						}

						for (int j=0; j<n_param; j++) {
							if (DS_integral.check(k,j+n_ode+n_ctrl)) {

								if (i>0)
									DS_obj(0, n_dis * (n_ode+n_ctrl) + j) += (T[i]-T[i-1])/2 * DS_integral(k,j+n_ode+n_ctrl)*lagrange_weight[k];
								if (i<n_dis-1)
									DS_obj(0, n_dis * (n_ode+n_ctrl) + j) += (T[i+1]-T[i])/2 * DS_integral(k,j+n_ode+n_ctrl)*lagrange_weight[k];
							}
						}
					}
				}

			}

			else if (twdiscretization->type==TWdiscretizationType::Trapez) {
				for (int i=0; i<n_dis; i++) {

					double *X1 = &X[(n_ode+n_ctrl)* i];
					double *U1 = &X[(n_ode+n_ctrl)* i + n_ode];

					double *P = &X[(n_ode+n_ctrl)* twdiscretization->stuetzstellen(n_dis)];

					phase->integral_diff(DS_integral,T[i],X1,U1,P);

					for (int k=0; k<n_integral; k++) {
						for (int j=0; j<n_ode + n_ctrl; j++) {
							if (DS_integral.check(k,j)) {

								if (i>0)
									DS_obj(0, i * (n_ode+n_ctrl) + j) += (T[i]-T[i-1])/2 * DS_integral(k,j)*lagrange_weight[k];
								if (i<n_dis-1)
									DS_obj(0, i * (n_ode+n_ctrl) + j) += (T[i+1]-T[i])/2 * DS_integral(k,j)*lagrange_weight[k];
							}
						}

						for (int j=0; j<n_param; j++) {
							if (DS_integral.check(k,j+n_ode+n_ctrl)) {

								if (i>0)
									DS_obj(0, n_dis * (n_ode+n_ctrl) + j) += (T[i]-T[i-1])/2 * DS_integral(k,j+n_ode+n_ctrl)*lagrange_weight[k];
								if (i<n_dis-1)
									DS_obj(0, n_dis * (n_ode+n_ctrl) + j) += (T[i+1]-T[i])/2 * DS_integral(k,j+n_ode+n_ctrl)*lagrange_weight[k];
							}
						}
					}
				}

			}


			else if (twdiscretization->type==TWdiscretizationType::HermiteSimpson) {
				for (int i=0; i<twdiscretization->stuetzstellen(n_dis); i++) {

					double *X1 = &X[(n_ode+n_ctrl)* i];
					double *U1 = &X[(n_ode+n_ctrl)* i + n_ode];

					double *P = &X[(n_ode+n_ctrl)* twdiscretization->stuetzstellen(n_dis)];

					if (i%2)
						phase->integral_diff(DS_integral,(T[(i-1)/2]+T[(i+1)/2])/2,X1,U1,P);
					else
						phase->integral_diff(DS_integral,T[i/2],X1,U1,P);


					for (int k=0; k<n_integral; k++) {
						for (int j=0; j<n_ode + n_ctrl; j++) {
							if (DS_integral.check(k,j)) {

								if (i%2) {
									DS_obj(0, i * (n_ode+n_ctrl) + j) += (T[(i+1)/2]-T[(i-1)/2])/6 * 4 * DS_integral(k,j)*lagrange_weight[k];
								} else {
									if (i/2>0)
										DS_obj(0, i * (n_ode+n_ctrl) + j) += (T[i/2]-T[i/2-1])/6 * DS_integral(k,j)*lagrange_weight[k];
									if (i/2<n_dis-1)
										DS_obj(0, i * (n_ode+n_ctrl) + j) += (T[i/2+1]-T[i/2])/6 * DS_integral(k,j)*lagrange_weight[k];
								}
							}
						}

						for (int j=0; j<n_param; j++) {
							if (DS_integral.check(k,j+n_ode+n_ctrl)) {

								if (i%2) {
									DS_obj(0, twdiscretization->stuetzstellen(n_dis) * (n_ode+n_ctrl) + j) += (T[(i+1)/2]-T[(i-1)/2])/6 * 4 * DS_integral(k,j+n_ode+n_ctrl)*lagrange_weight[k];
								}else {
									if (i/2>0)
									DS_obj(0, twdiscretization->stuetzstellen(n_dis) * (n_ode+n_ctrl) + j) += (T[i/2]-T[i/2-1])/6 * DS_integral(k,j+n_ode+n_ctrl)*lagrange_weight[k];
								if (i/2<n_dis-1)
									DS_obj(0, twdiscretization->stuetzstellen(n_dis) * (n_ode+n_ctrl) + j) += (T[i/2+1]-T[i/2])/6 * DS_integral(k,j+n_ode+n_ctrl)*lagrange_weight[k];

								}

							}
						}

					}
				}

				/**for (int i=0;i<n_param;i++) {
					if (DS_integral.check(k,n_ode + n_ctrl + i)) {

						//DS_obj(0, n_dis * (n_ode+n_ctrl) + i); ??????? TODO

					}
				}*/
			}

			else if (twdiscretization->type==TWdiscretizationType::Lobatto) {

				cout << "TODO 3 Lobatto" << endl;
			/*	for (int i=0; i<twdiscretization->stuetzstellen(n_dis); i++) {

					double *X1 = &X[(n_ode+n_ctrl)* i];
					double *U1 = &X[(n_ode+n_ctrl)* i + n_ode];

					double *P = &X[(n_ode+n_ctrl)* twdiscretization->stuetzstellen(n_dis)];

					if (i%2)
						integral_diff(DS_integral,(T[(i-1)/2]+T[(i+1)/2])/2,X1,U1,P);
					else
						integral_diff(DS_integral,T[i/2],X1,U1,P);


					for (int k=0; k<n_integral; k++) {
						for (int j=0; j<n_ode + n_ctrl; j++) {
							if (DS_integral.check(k,j)) {

								if (i%2) {
									DS_obj(0, i * (n_ode+n_ctrl) + j) += (T[(i+1)/2]-T[(i-1)/2])/6 * 4 * DS_integral(k,j)*lagrange_weight[k];
								} else {
									if (i/2>0)
										DS_obj(0, i * (n_ode+n_ctrl) + j) += (T[i/2]-T[i/2-1])/6 * DS_integral(k,j)*lagrange_weight[k];
									if (i/2<n_dis-1)
										DS_obj(0, i * (n_ode+n_ctrl) + j) += (T[i/2+1]-T[i/2])/6 * DS_integral(k,j)*lagrange_weight[k];
								}
							}
						}

						for (int j=0; j<n_param; j++) {
							if (DS_integral.check(k,j+n_ode+n_ctrl)) {

								if (i%2) {
									DS_obj(0, twdiscretization->stuetzstellen(n_dis) * (n_ode+n_ctrl) + j) += (T[(i+1)/2]-T[(i-1)/2])/6 * 4 * DS_integral(k,j+n_ode+n_ctrl)*lagrange_weight[k];
								}else {
									if (i/2>0)
									DS_obj(0, twdiscretization->stuetzstellen(n_dis) * (n_ode+n_ctrl) + j) += (T[i/2]-T[i/2-1])/6 * DS_integral(k,j+n_ode+n_ctrl)*lagrange_weight[k];
								if (i/2<n_dis-1)
									DS_obj(0, twdiscretization->stuetzstellen(n_dis) * (n_ode+n_ctrl) + j) += (T[i/2+1]-T[i/2])/6 * DS_integral(k,j+n_ode+n_ctrl)*lagrange_weight[k];

								}

							}
						}

					}
				}*/

				/**for (int i=0;i<n_param;i++) {
					if (DS_integral.check(k,n_ode + n_ctrl + i)) {

						//DS_obj(0, n_dis * (n_ode+n_ctrl) + i); ??????? TODO

					}
				}*/
			}

		}

	}



	if (ret) {
		for (int k=0; k<DF_nnz; k++) {
			DF.val[DF_start + k] = ScaleObj * DS_obj.get(0,DF.row[DF_start+k]-1-Delta1);
		}
		twcount_calls.obj.df++;
		return;
	}

	// Ableitung numerisch
	for (int k=0; k<DF_nnz; k++) {

		int index = DF.row[DF_start + k]-1 - Delta1;

		//cout << Delta1 << "   " << X << " " << index << " " << k << " " <<  DF_start << endl;

		double v = X[index];
		X[index] += eps;

		twcount_calls.obj.f+=2;
		Lagrange();
		double f1 = phase->obj();
		for (int i=0; i<n_integral; i++) {
			f1 += lagrange_weight[i] * lagrange_integral[n_integral*(twdiscretization->stuetzstellen(n_dis)-1)+i];
		}

		X[index] = v- eps;
		Lagrange();
		double f2 = phase->obj();
		for (int i=0; i<n_integral; i++) {
			f2 += lagrange_weight[i] * lagrange_integral[n_integral*(twdiscretization->stuetzstellen(n_dis)-1)+i];
		}

		// zentraler Differenzenquotient
		DF.val[DF_start + k] = ScaleObj * (f1-f2)/(2*eps);

		X[index] = v;
	}
}




std::ostream& operator<<(std::ostream &os, const FullDisTransWorhp &p) {

	os << "Phase [" << p.n_ode << " " << p.n_ctrl << " " << p.n_dis << " "
	   << p.n_param << "] = " << p.n_var <<  endl;

	if (p.X) {

		for (int i=0; i<p.n_ode; i++) {
			os << std::setw(13) << "X" << i;
		}
		for (int i=0; i<p.n_ctrl; i++) {
			os << std::setw(13) << "U" << i;
		}
		os << endl;

		for (int j=0; j<p.n_dis; j++) {

			for (int i=0; i<p.n_ode; i++) {
				os << std::setw(14) << p.x(j,i);
			}

			for (int i=0; i<p.n_ctrl; i++) {
				os << std::setw(14) << p.u(j,i);
			}
			os << endl;
		}

	}
// if (p.next) os << *p.next;

	return os;
}





void FullDisTransWorhp::init0() {

	n_var = (n_ode + n_ctrl) * twdiscretization->stuetzstellen(n_dis) + n_param;
	n_con = n_ode * twdiscretization->stufen * (n_dis-1) + n_rand + n_neben*twdiscretization->stuetzstellen(n_dis);
	if (n_con<0) n_con = 0;

	tmp_gg_1 = vector<double>(n_con);
	tmp_gg_2 = vector<double>(n_con);
	tmp_gg_3 = vector<double>(n_con);
	tmp_gg_4 = vector<double>(n_con);
	tmp_gg_5 = vector<double>(n_con);

	phase->localinit();

	DS_obj.Init(1,         n_var);
	DS_ode.Init(n_ode,     n_ode+n_ctrl + n_param);
	DS_rand.Init(n_rand,   n_var);
	DS_neben.Init(n_neben, n_ode+n_ctrl + n_param);

	DS_integral.Init(n_integral, n_ode+n_ctrl + n_param);
}


void FullDisTransWorhp::Debug_G() {

	/*for (int i=0;i<n_con;i++) {
	 cout << setw(14) << G_low[i] << setw(14) << G[i] << setw(14) << G_upp[i];

	 if (i<n_ode*n_dis) {
	  cout << "  ODE:  " << setw(4) << (i%n_dis) << setw(4) << (i/n_dis);
	 }
	 else if (i<n_ode*n_dis + n_rand) {
	  cout << "  RAND: " << (i-n_ode*n_dis);
	 }
	 else {
	  int ii= (i-n_ode*n_dis-n_rand);
	  cout << "  NEBEN:" << (ii%n_dis) << " " << (ii/n_dis);

	 }
	 cout<< endl;
	}*/

	double mm = 0;
	int mi = 0;
	int mj = 0;

	for (int i=0; i<n_dis-1; i++) {

		cout << "ODE: " << std::setw(4) << (i);
		for (int j=0; j<n_ode; j++) {
			cout << std::setw(14) << G[i*n_ode+j];

			if (std::abs(G[i*n_ode+j]) > mm) {
				mm=std::abs(G[i*n_ode+j]);
				mi=i;
				mj=j;
			}
		}
		cout << endl;
	}

	cout << "Max"  << std::setw(14)<< mm << std::setw(4) << mi << std::setw(4) << mj << endl;

}



void FullDisTransWorhp::Debug_G2() {

	double mm = 0;
	int mi = 0;
	int mj = 0;
	for (int i=0; i<n_dis-1; i++) {

		for (int j=0; j<n_ode; j++) {

			if (std::abs(G[i*n_ode+j]) > mm) {
				mm=std::abs(G[i*n_ode+j]);
				mi=i;
				mj=j;
			}
		}
	}

	cout << "Max"  << std::setw(14)<< mm << std::setw(4) << mi << std::setw(4) << mj << endl;
}

std::string FullDisTransWorhp::type_G(int row) const {

	stringstream a;

	if (row < n_ode*(twdiscretization->stuetzstellen(n_dis)-1)) {
		a << " ODE  " << std::setw(2) << row%n_ode << " ";
	}
	else if (row < n_ode*(twdiscretization->stuetzstellen(n_dis)-1)+n_rand) {
		a << " RAND " << std::setw(2) << (row-n_ode*(twdiscretization->stuetzstellen(n_dis)-1))%n_rand << " ";
	}
	else if (row < n_ode*(twdiscretization->stuetzstellen(n_dis)-1)+n_rand+n_neben*(twdiscretization->stuetzstellen(n_dis))) {
		a << " NB   " << std::setw(2) << (row-n_ode*(twdiscretization->stuetzstellen(n_dis)-1)-n_rand-n_neben)%n_neben << " ";
	}

	return a.str();
}


void FullDisTransWorhp::fromMATLAB_impl(std::ifstream &stream) {

	std::string line;

	/// Header
	std::getline(stream, line);
//cout << ":" << line << endl;
	/// Parameter p
	std::getline(stream, line);
	std::string line2 = std::string(line,2);
	std::vector<double> v = ToDoubleArray(line2);
	for (int i=0; i<n_param; i++) {
		X[p_index(i)] = v[i];
	}

	/// Data
	std::getline(stream, line);
//cout << ":" << line << endl;
	std::vector<double> v_last;

	v = ToDoubleArray(line);
	if (phase->freetime) v[0] = v[0]/p(0);

	v_last = v;
//cout << ":" << v_last << " " << n_dis << endl;
	for (int i=0; i<n_dis; i++) {

		while (T[i] > v[0]) {
//cout << i << "   " << T[i] << " " << v[0] << endl;
			std::getline(stream, line);

			if (!stream) {
				v_last = v;
				break;
			} else {
				v_last = v;

				v = ToDoubleArray(line);
				if (phase->freetime) v[0] = v[0]/p(0);
			}
		}

		double vv=1;

		if (v[0]!=v_last[0])
			vv = (T[i] - v_last[0]) / (v[0] - v_last[0]);
//cout << i << "   vv " << vv << endl;

		for (int j=0; j<n_ode; j++) {
			X[x_index(i,j)] = v_last[j+1] + (v[j+1]-v_last[j+1]) * vv;

			if (i>0 && twdiscretization->type==TWdiscretizationType::HermiteSimpson) {
				X[x_index__(2*i-1,j)] = (X[x_index__(2*(i-1),j)] + X[x_index__(2*i,j)])/2;
			}
		}

		for (int j=0; j<n_ctrl; j++) {
			X[u_index(i,j)] = v_last[j+n_ode+1] + (v[j+n_ode+1]-v_last[j+n_ode+1]) * vv;

			if (i>0 && twdiscretization->type==TWdiscretizationType::HermiteSimpson) {
				X[u_index__(2*i-1,j)] = (X[u_index__(2*(i-1),j)] + X[u_index__(2*i,j)])/2;
			}
		}

	}

#ifdef TRANSWORHP_GRAPHICS
	if (viewer) viewer->autoScale();
#endif
}


void FullDisTransWorhp::GetState(double *xx, double t) {

	int i=0;
	for (; i<n_dis; i++) {
		if (t < T[i]) break;
	}

	for (int j=0; j<n_ode; j++) {
		if (i>0 && i<n_dis)
			xx[j] = x(i-1,j) +  (x(i,j) - x(i-1,j)) / (T[i]-T[i-1])*(t-T[i-1]);
		if (i>=n_dis)
			xx[j] = x(n_dis-1,j);
		if (i==0)
			xx[j] = x(0,j);
	}

	//cout << t << " " << 0 << " " << xx[0] << endl;
}


void FullDisTransWorhp::GetControl(double *uu, double t) {

	int i=0;
	for (; i<n_dis; i++) {
		if (t < T[i]) break;
	}

	for (int j=0; j<n_ctrl; j++) {
		if (i>0 && i<n_dis)
			uu[j] = u(i-1,j) +  (u(i,j) - u(i-1,j)) / (T[i]-T[i-1])*(t-T[i-1]);
		if (i>=n_dis)
			uu[j] = u(n_dis-1,j);
		if (i==0)
			uu[j] = u(0,j);
	}

	//cout << t << " " << 0 << " " << xx[0] << endl;
}


#ifdef TRANSWORHP_GRAPHICS
string FullDisTransWorhp::setViewerTitle() {
	return string(" [ndis = " + std::to_string(n_dis) + "]");
}


void FullDisTransWorhp::updateViewer(DataStorage* ds, vector<double> &tmptime) {

	ds->setData( 	twdiscretization->stuetzstellen(n_dis),tmptime.data(),
			X, n_ode+n_ctrl, n_ode,
			&X[n_ode], n_ode+n_ctrl, n_ctrl,
			G,  0,  0,
			lagrange_integral.data(),  n_integral,n_integral);

	if (twparameter->showGrid) {
		if (twdiscretization->type == TWdiscretizationType::Euler || twdiscretization->type == TWdiscretizationType::Trapez) {
			ds->addDisData(n_dis, n_dis, T.data(), X, &X[n_ode]);
		} else if (twdiscretization->type == TWdiscretizationType::HermiteSimpson) {
			ds->addDisData(n_dis, n_dis, 2*(n_ode+n_ctrl), 2*(n_ode+n_ctrl), T.data(), X, &X[n_ode]);
		}
	}

}


void FullDisTransWorhp::setTemptimeForViewer(vector<double> &tmptime) {

	if (twdiscretization->type==TWdiscretizationType::Trapez) {
		tmptime = T;
	}
	else if (twdiscretization->type==TWdiscretizationType::HermiteSimpson) {
		tmptime.resize(n_dis*2);
		for (int i=0;i<n_dis;i++) {
			tmptime[2*i] = T[i];
			if (i<n_dis-1) tmptime[2*i+1] = (T[i] + T[i+1])/2.;
		}
	}
	else if (twdiscretization->type==TWdiscretizationType::Lobatto) {
		tmptime.resize(n_dis*3);
		for (int i=0;i<n_dis;i++) {
			tmptime[3*i] = T[i];
			if (i<n_dis-1) {
				tmptime[3*i+1] = T[i] + (T[i+1]-T[i]) * (5.-sqrt(5.))/10.;
				tmptime[3*i+2] = T[i] + (T[i+1]-T[i]) * (5.+sqrt(5.))/10.;
			}
		}
	}
	else if (twdiscretization->type==TWdiscretizationType::Euler) {
		tmptime = T;
	}
	else if (twdiscretization->type==TWdiscretizationType::MultipleShooting) {
		tmptime = T;
	}
}
#endif


int FullDisTransWorhp::Integrate(int btableau) {

	int steps = 0;

	if (btableau==6) {

		if (butcher.stufen()==0) butcher.Init(btableau,twparameter->stepsize,false);

		int startflag = 1;
		for (int i=0; i<n_dis-1; i++) {
			int a = integrateRKF(i, startflag);
			startflag=2;
			steps+=a;

		//	cout << setw(5) << i << " " << setw(10) << a << endl;
		}

	}
	else {

		if (butcher.stufen()==0) butcher.Init(btableau,twparameter->stepsize);

		for (int i=0; i<n_dis-1; i++) {
			int a = integrate(i);
			steps+=a;

		//	cout << setw(5) << i << " " << setw(10) << a << endl;
		}

	}

	if (twdiscretization->type==TWdiscretizationType::HermiteSimpson) {

		// Zwischenwerte interpolieren fuer schoenere Optik...
		// allgemein immer machen wenn #Zwischenpunkte>1
		for (int i=1; i<n_dis; i++) {
			for (int l=0; l<n_ode; l++) {
				X[x_index(i-1,l) + n_ode+n_ctrl] = .5*(X[x_index(i-1,l)] + X[x_index(i,l)]);
			}
			for (int l=0; l<n_ctrl; l++) {
				X[u_index(i-1,l) + n_ode+n_ctrl] = .5*(X[u_index(i-1,l)] + X[u_index(i,l)]);
			}
		}
	}

	else if (twdiscretization->type==TWdiscretizationType::Lobatto) {

		// Zwischenwerte interpolieren fuer schoenere Optik...
		// allgemein immer machen wenn #Zwischenpunkte>1
		for (int i=1; i<n_dis; i++) {
			for (int l=0; l<n_ode; l++) {
				X[x_index(i-1,l) + n_ode+n_ctrl] = X[x_index(i-1,l)] + (X[x_index(i,l)]- X[x_index(i,l-1)]) * (5-sqrt(5.)/10);
			}
			for (int l=0; l<n_ctrl; l++) {
				X[u_index(i-1,l) + n_ode+n_ctrl] = X[x_index(i-1,l)] + (X[x_index(i,l)]- X[x_index(i,l-1)]) * (5+sqrt(5.)/10);
			}
		}
	}

#ifdef TRANSWORHP_GRAPHICS
	if (viewer) viewer->autoScale();
#endif

	return steps;

}

int FullDisTransWorhp::integrateRKF(int index, int &startflag) {

	double t  = T[index];
	double tend = T[index+1];

	for (int l=0; l<n_ode; l++) {
		tmp_ode_1[l] = x(index,l);
	}
	for (int l=0; l<n_ctrl; l++) {
		tmp_ctrl_1[l] = u(index,l);
		tmp_ctrl_2[l] = u(index+1,l);
	}

	int a = butcher.RungeKuttaF(phase,t,T[index],tend,tmp_ode_1.data(),tmp_ctrl_1.data(),tmp_ctrl_2.data(),&X[p_index(0)],startflag);

	for (int l=0; l<n_ode; l++){
		setx(index+1,l) = tmp_ode_1[l];
	}

	return a;

}

int FullDisTransWorhp::integrate(int index) {

	double h  = .001;
	double t  = T[index];
	double tend = T[index+1];

	int ct = 0;

	for (int l=0; l<n_ode; l++) {
		tmp_ode_1[l] = x(index,l);
	}
	for (int l=0; l<n_ctrl; l++) {
		tmp_ctrl_1[l] = u(index,l);
		tmp_ctrl_2[l] = u(index+1,l);
	}

	while (t<tend) {

	//	cout << "t<tend" << t << " " << tend << endl;
		if (h+t>tend) {
			h = tend-t;
		}

		//cout << t << " " << h << endl;
		ct++;

		int a = butcher.RungeKutta(phase,t,T[index],tend,tmp_ode_1.data(),tmp_ctrl_1.data(),tmp_ctrl_2.data(),&X[p_index(0)],h);


		if (a==-1) {
			cout << "Too small at " << t  << endl;
			break;
		}

		if (a==0) {
			if (t==tend) {

				for (int l=0; l<n_ode; l++) {
					setx(index+1,l) = tmp_ode_1[l];
				}
			}
		}
	}

	return ct;
}


void FullDisTransWorhp::GetBoundaryIndices(vector<int> &indices, int d) {

	if ( X_low[x_index__(0,d)] == X_upp[x_index__(0,d)] ) {
		int ind = twdiscretization->stuetzstellen(1) ;
		indices.push_back(ind-1);
		//cout << it-twfolder->phases.begin() << " " << s<< " L " << ii << " " << ind-1 << endl;
	}
	if ( X_low[x_index(n_dis-1,d)] == X_upp[x_index(n_dis-1,d)] ) {
		int ind = twdiscretization->stuetzstellen(n_dis) ;//+ ct;
		indices.push_back(ind-1);
		//cout << it-twfolder->phases.begin() << " " << s << " R " << ii << " " << ind-1 << endl;
	}
}

/*
Matthias Rick
Hilfsfunktion
*/
double FullDisTransWorhp::auxFehlerInt(double t, int komp, TWbaseSpline **zustand, TWbaseSpline **steuerung) {

	// nutze existierende Pointer
	auto DX = tmp_ode_2.data(); //ode-Wert
	auto state = tmp_ode_1.data(); //Zustaende
	auto ctrl = tmp_ctrl_1.data(); //Steuerung

	for (int i = 0; i < n_ode; i++) {
		state[i] = zustand[i]->eval(t);
	}

	for (int i = 0; i < n_ctrl; i++) {
		ctrl[i] = steuerung[i]->eval(t);
	}

	auto P = &X[(n_ode + n_ctrl)* twdiscretization->stuetzstellen(n_dis)]; //Parameter

	phase->ode(DX, t, state, ctrl, P);

	return std::abs(zustand[komp]->deriv(t) - DX[komp]);
}



/*
Matthias Rick
bestimmt den Diskretisierungsfehler im Intervall [a,b] fuer alle Zustaende
berechnet Integral ( |spline(x)' - f(spline(x),spline(u),t)| )
*/
void FullDisTransWorhp::disFehlerIntegral(double *fehler, const double a, const double b, TWbaseSpline **zustand, TWbaseSpline **steuerung) {
	const int tiefe = 7; // Genauigkeit des Verfahren
	double Tra[tiefe][tiefe]; // eigentlich Dreieckmatrix
	double h[tiefe]; // Schrittweite

	//Schrittweiten berechnen
	h[0] = (b - a);
	for (int i = 1; i < tiefe; i++) {
		h[i] = h[0] / pow(2., i);
	}

	for (int komp = 0; komp < n_ode; komp++) {

		//ersten Wert per Trapezverfahren
		Tra[0][0] = 0.5*h[0] * (auxFehlerInt(a, komp, zustand, steuerung) + auxFehlerInt(b, komp, zustand, steuerung));

		//restliche Werte nach zusammengesetzter Trapezregel
		for (int i = 1; i < tiefe; i++) {

			// Schleife fuer Summe
			double sum = 0;
			for (int j = 0; j <= pow(2., i - 1) - 1; j++) {
				sum += auxFehlerInt(a + 0.5*h[i - 1] + j * h[i - 1], komp, zustand, steuerung);
			}

			Tra[0][i] = 0.5*Tra[0][i - 1] + 0.5*h[i - 1] * sum;
		}

		//weitere Spalten nach Rombergverfahren (Bronstein S.972)
		for (int i = 1; i < tiefe; i++) {
			for (int k = 1; k <= i; k++) {
				Tra[k][i] = Tra[k - 1][i] + (Tra[k - 1][i] - Tra[k - 1][i - 1]) / (pow(4., k) - 1);
			}
		}

		// Ausgabe des Fehlers
		//cout << "Intervall [" << a << "," << b << "]" << "komp=" << komp << " : Fehler=" << Tra[tiefe - 1][tiefe - 1] << endl;

		fehler[komp] = Tra[tiefe - 1][tiefe - 1];


		//ausgabe des schemas
		/*
		for (int i = 0; i < tiefe; i++) {
			for (int k = 0; k <= i; k++) {
				cout << Tra[k][i] << " \t";
			}
			cout << endl;
		}
		*/
	}
}

/*
Matthias Rick
berechnet den Diskretisierungsfehler nach Betts
*/
void FullDisTransWorhp::diskretisierungsfehlerBetts() {

	auto DX1 = new double[n_ode]; // Wert rechte Seite (linkte Stuetzstelle)
	auto DX2 = new double[n_ode]; // Wert rechte Seite (rechte Stuetzstelle)
	auto DX05 = new double[n_ode]; // Wert rechte Seite (mittige Stuetzstelle)

	auto fehler = new double[n_ode]; // zur Uebergabe an disFehlerIntegral()

	auto eta = new double*[n_ode]; // eta_i,k (Betts 4.154)
	for (int i = 0; i < n_ode; i++) {
		eta[i] = new double[n_dis];
	}

	auto hilf = new double*[n_ode]; // Hilfsvariable fuer Bestimmung der Gewichte (Ableitung von Spline wird drin gespeichert)
	for (int i = 0; i < n_ode; i++) {
		hilf[i] = new double[n_dis];
	}

	double *X05=nullptr;
	double *U05=nullptr;

	// Splines erzeugen
	auto zustand = new TWbaseSpline*[n_ode];
	auto steuerung = new TWbaseSpline*[n_ctrl];

	for (int dis = 0; dis < n_dis - 1; dis++) {

		double xk = T[dis];
		double xk1 = T[dis+1];

		auto X1 = &X[(n_ode + n_ctrl)* twdiscretization->punkte()*dis];
		auto U1 = &X[(n_ode + n_ctrl)* twdiscretization->punkte()*dis + n_ode];

		auto X2 = &X[(n_ode + n_ctrl)* twdiscretization->punkte()*(dis+1)];
		auto U2 = &X[(n_ode + n_ctrl)* twdiscretization->punkte()*(dis+1) + n_ode];

		if (twdiscretization->type == TWdiscretizationType::HermiteSimpson) {
			X05 = &X[(n_ode+n_ctrl)* (twdiscretization->punkte()*dis+1)];
			U05 = &X[(n_ode+n_ctrl)* (twdiscretization->punkte()*dis+1) + n_ode];
		}

		auto P = &X[(n_ode + n_ctrl)* twdiscretization->stuetzstellen(n_dis)];

		// rechte Seite auswerten
		phase->ode(DX1, xk, X1, U1, P);
		phase->ode(DX2, xk1,X2, U2, P);
		if (twdiscretization->type == TWdiscretizationType::HermiteSimpson) {
			phase->ode(DX05, (xk+xk1)/2., X05, U05, P);
		}

		for (int i = 0; i < n_ode; i++) {

			if (twdiscretization->type != TWdiscretizationType::HermiteSimpson) {
				zustand[i] = new TWspline(xk, xk1, x(dis, i), x(dis + 1, i), DX1[i], DX2[i]);
			} else {
				zustand[i] = new TWsplineHS(xk, (xk+xk1)/2., xk1, x(dis, i), x__(2*dis+1, i), x(dis+1, i), DX1[i], DX05[i], DX2[i]);
			}

			//fuer Gewichtung - Wert der Ableitung
			hilf[i][dis] = std::abs(zustand[i]->deriv(xk));
		}

		for (int i = 0; i < n_ctrl; i++) {

			if (twdiscretization->type != TWdiscretizationType::HermiteSimpson) {
				steuerung[i] = new TWspline(xk, xk1, u(dis, i), u(dis+1, i));
			} else {
				steuerung[i] = new TWsplineHS(xk, (xk+xk1)/2., xk1, u(dis, i), u__(2*dis+1, i), u(dis+1, i));
			}
		}

		// Fehler berechnen
		disFehlerIntegral(fehler, xk, xk1, zustand, steuerung);

		// eta berechnen
		for (int i = 0; i < n_ode; i++) {
			eta[i][dis] = fehler[i];
		}

		//cout << "...................................." << endl;

		// Speicher freigeben
		for (int i = 0; i < n_ode; i++) {
			delete zustand[i];
		}
		for (int i = 0; i < n_ctrl; i++) {
			delete steuerung[i];
		}

	}

	// Speicher freigeben
	delete[] zustand;
	delete[] steuerung;

	// Gewichtung w_i bestimmen Betts (4.157)
	//cout << "Gewichtung:" << endl;
	auto gewichtung = new double[n_ode]; // Gewichtung der Fehler
	for (int i = 0; i < n_ode; i++) {
		gewichtung[i] = 0;
		for (int k = 0; k < n_dis - 1; k++) {
			if (std::max(hilf[i][k], x(i, k)) > gewichtung[i]) {
				gewichtung[i] = std::max(hilf[i][k], x(i, k));
			}
		}
		//cout << "w_" << i << " = " << gewichtung[i] << endl;
	}

	//cout << "---------------------................------------------" << endl;

	auto disFehler = new double[n_dis]; // Fehler eps_k (Betts 4.156)

	vector<double> fehler1; // Rueckgabewert
	fehler1.reserve(n_dis); //Groesse festlegen

	// realativer lokaler Fehler Betts (4.156)
	for (int k = 0; k < n_dis - 1; k++) {
		disFehler[k] = 0;
		for (int i = 0; i < n_ode; i++) {
			double aux = eta[i][k] / (gewichtung[i] + 1);
			if (aux > disFehler[k]) {
				disFehler[k] = aux;
			}
		}
		//cout << "eps_" << k << " = " << disFehler[k] << endl;
		fehler1.push_back(disFehler[k]);
	}

	stringstream sstm; // Textausgabe

	MyStatus("MeshRef", "error mode: Betts", Status::NORMAL);

	// durchschnittlicher Fehler Betts (4.169)
	double aveFehler = 0;
	for (int i = 0; i < n_dis - 1; i++) {
		aveFehler += disFehler[i];
	}
	aveFehler = aveFehler / (n_dis - 1);
	//cout << "durschnittlicher Fehler: " << aveFehler << endl;
	sstm << "average error: " << aveFehler;
	MyStatus("MeshRef", sstm.str(), Status::NORMAL);
	sstm.str("");

	// Intervall mit dem groessten Fehler finden
	double max = disFehler[0];
	int intAux = 0;
	for (int i = 1; i < n_dis-1; i++) {
		if (disFehler[i] > max) {
			max = disFehler[i];
			intAux = i;
		}
	}
	//cout << "Intervall mit groesstem Fehler: " << intAux << " mit: " << max << endl;
	sstm << "interval with largest error: " << intAux << " with: " << max;
	MyStatus("MeshRef",sstm.str(), Status::NORMAL);
	sstm.str("");

	// Speicher freigeben
	delete[] DX1;
	delete[] DX2;
	delete[] DX05;
	delete[] fehler;
	delete[] disFehler;
	delete[] gewichtung;

	for (int i = 0; i < n_ode; i++) {
		delete[] eta[i];
		delete[] hilf[i];
	}
	delete[] eta;
	delete[] hilf;

	FEHLER.push_back(std::move(fehler1));
}

/*
Matthias Rick
berechnet den Diskretisierungsfehler ueber verschiedene (hoehere Ordnung) Diskretisierungen
*/
void FullDisTransWorhp::diskretisierungsfehler2() {

	vector<double> fehler;
	fehler.reserve(n_dis);

	double *X1;
	double *X2;
	double *U1;
	double *U2;
	double *P = &X[(n_ode+n_ctrl)* twdiscretization->stuetzstellen(n_dis)];

	// euler -> trapez
	if (twdiscretization->type==TWdiscretizationType::Euler) {

		auto x1 = new double[n_ode];

		for (int i = 0; i < n_dis-1; i++) {

			X1 = &X[(n_ode+n_ctrl)* i];
			U1 = &X[(n_ode+n_ctrl)* i + n_ode];

			phase->ode(tmp_ode_1.data(),T[i],X1,U1,P);

			X2 = &X[(n_ode+n_ctrl)* (i+1)];
			U2 = &X[(n_ode+n_ctrl)* (i+1) + n_ode];

			phase->ode(tmp_ode_2.data(),T[i+1],X2,U2,P);

			double err = 0;

			for (int k = 0; k < n_ode; k++) {
				x1[k] =  X[x_index(i,k)] + (T[i+1]-T[i])/2.0 * (tmp_ode_1[k] + tmp_ode_2[k]);
				if (std::abs(x1[k] - X[x_index(i+1,k)]) > err) {
					err = std::abs(x1[k] - X[x_index(i+1,k)]);
				}
			}

			fehler.push_back(err);
		}

		delete[] x1;
	}

	// trapez -> hermite simpson
	else if (twdiscretization->type == TWdiscretizationType::Trapez) {

		auto X05 = new double[n_ode];
		auto U05 = new double[n_ctrl];
		auto x1 = new double[n_ode];

		auto DX05 = new double[n_ode];

		for (int i = 0; i < n_dis-1; i++) {

			double err = 0;

			X1 = &X[(n_ode+n_ctrl)* i];
			U1 = &X[(n_ode+n_ctrl)* i + n_ode];

			phase->ode(tmp_ode_1.data(),T[i],X1,U1,P);

			X2 = &X[(n_ode+n_ctrl)* (i+1)];
			U2 = &X[(n_ode+n_ctrl)* (i+1) + n_ode];

			phase->ode(tmp_ode_2.data(),T[i+1],X2,U2,P);

			//cout << "x_" << i <<"+1/2 berechnen" << endl;
			// x_i+1/2 berechnen
			for (int k = 0; k < n_ode; k++) {
				X05[k] = 0.5*(X[x_index(i+1,k)] + X[x_index(i,k)]) + (T[i+1]-T[i])/8.0*(tmp_ode_1[k] - tmp_ode_2[k]);
				//cout << x05[k] << " ";
			}
			//cout << endl;

			// u_i+1/2 berechnen
			for (int k = 0; k < n_ctrl; k++) {
				U05[k] = 0.5*(X[u_index(i,k)]+X[u_index(i+1,k)]);
			}

			phase->ode(DX05,0.5*(T[i]+T[i+1]),X05,U05,P);

			//cout << "x_" << i << "+1 berechnen" << endl;

			//x_i+1 berechnen
			for (int k = 0; k < n_ode; k++) {
				x1[k] =  X[x_index(i,k)] + (T[i+1]-T[i])/6.0*(tmp_ode_2[k] + 4.0*DX05[k] + tmp_ode_1[k]);
				if (std::abs(x1[k] - X[x_index(i+1,k)]) > err) {
					err = std::abs(x1[k] - X[x_index(i+1,k)]);
				}

				//cerr << std::abs(x1[k] - X[x_index(i+1,k)]) << "\t";
			}
			fehler.push_back(err);
			//cerr << endl;
		}

		//cout << "Fehler mit hoeheren Verfahren" << endl;
		MyStatus("MeshRef", "error mode: high order", Status::NORMAL);

		// durchschnittlicher Fehler
		double aveFehler = 0;
		for (int i = 0; i < n_dis - 1; i++) {
			aveFehler += fehler[i];
		}
		aveFehler = aveFehler / (n_dis - 1);
		//cout << "durschnittlicher Fehler: " << aveFehler << endl;


		string s( "average error: " + std::to_string(aveFehler) );
		MyStatus("MeshRef", s, Status::NORMAL);

		// Intervall mit dem groessten Fehler finden
		double max = fehler[0];
		int intAux = 0;
		for (int i = 1; i < n_dis-1; i++) {
			if (fehler[i] > max) {
				max = fehler[i];
				intAux = i;
			}
		}
		//cout << "Intervall mit groesstem Fehler: " << intAux << " mit: " << max << endl;
		s = "interval with largest error: " + std::to_string(intAux) + " with: " + std::to_string(max);
		MyStatus("MeshRef", s, Status::NORMAL);

		delete[] X05;
		delete[] U05;
		delete[] x1;

		delete[] DX05;
	}

	FEHLER.push_back(std::move(fehler));
}

/*
Matthias Rick
berechnet den Diskretisierungsfehler
*/
void FullDisTransWorhp::diskretisierungsfehler() {

	if (twdiscretization->type == TWdiscretizationType::HermiteSimpson) {
		diskretisierungsfehlerBetts();
	} else {
		if (twparameter->meshref_err_mod == 0) {
			diskretisierungsfehlerBetts();
		} else {
			diskretisierungsfehler2();
		}
	}
}

/*
Matthias Rick
berechnet neue Diskretisierung
*/
vector<double> FullDisTransWorhp::refineAlg(double &max) {

	// I initialisieren
	vector<int> I(n_dis-1);

	if (twparameter->meshref_mod == 0 || twparameter->meshref_mod == 2) {
		betts(I, max);
	}

	//Anzahl neuer Punkte bestimmen
	const int anz = std::accumulate(I.begin(),I.end(),0);

	vector<double> T_neu; // neues Gitter
	T_neu.reserve(anz);

	T_neu.push_back(T[0]);

	for (int i = 1; i < n_dis; i++) {

		for (int m = 1; m <= I[i - 1]; m++) {
			// Zwischenwerte interpolieren
			T_neu.push_back(static_cast<double>(m) / (I[i - 1] + 1) * (T[i] - T[i - 1]) + T[i - 1]);
		}

		T_neu.push_back(T[i]);
	}

	if (twparameter->meshref_mod == 1 || twparameter->meshref_mod == 2) {

		vector<double> bueskes( bueskens() ); // Vektor fuer Auf-/Absprungpunkte

		if (!bueskes.empty()) {
			T_neu.insert(T_neu.end(),bueskes.begin(),bueskes.end());
			std::sort(T_neu.begin(), T_neu.end());
		}

		// Punkte loeschen, die zu dicht liegen
		//cerr << "size: " << T_neu.size() << endl;
		for (size_t i = 1; i < T_neu.size(); i++) {
			//cerr << T_neu[i] - T_neu[i-1] << endl;
			//cerr << "i:" << i << " " << T_neu[i] << endl;
			if (T_neu[i] - T_neu[i-1] < 1e-6) {
				std::cerr << "loesche " << T_neu[i] << " da zu dicht an " << T_neu[i-1] << endl;

				// anderer Ansatz: Mittelwert der Punkte
				T_neu[i-1] = 0.5*(T_neu[i]+T_neu[i-1]);

				T_neu.erase(T_neu.begin()+i);
				i--;
			}
		}
	}

	return T_neu;
}

void FullDisTransWorhp::betts(vector<int> &I, double &max) {

	int M; // maximale Anzahl der Punkte die eingefuegt werden (so viele, wie Startpunkte)
	const int M1 = twparameter->meshref_M1; // maximale Anzahl der Punkte pro Intervall

	const double tol = twparameter->meshref_tol; // Toleranz

	int anzI = 0; // Summe der Eintraege von I

	const double kappa = twparameter->meshref_kappa;

	int ordVerf = 1; // Ordnung des Verfahrens
	int redu = twparameter->meshref_R; // Reduktion der Ordnung

	bool flag = false; // fuer Abbruchkriterium

	vector<double> fehlerAux; //zum zwischenspeichern des Fehlervektors

	// M setzen
	if (twparameter->meshref_M == -1) {
		M = n_dis - 1;
	} else {
		M = twparameter->meshref_M;
	}

	// Ordnung des Verfahrens
	if (twparameter->twdiscretization.type == TWdiscretizationType::Euler) {
		ordVerf = 1;
	} else if (twparameter->twdiscretization.type == TWdiscretizationType::Trapez) {
		ordVerf = 2;
	} else if (twparameter->twdiscretization.type == TWdiscretizationType::HermiteSimpson) {
		ordVerf = 4;
	}

	if (redu > ordVerf) {
		redu = ordVerf;
	}

	if (!FEHLER.empty()){
		fehlerAux = FEHLER.back();
	}

	/* Entfernen von Punkten
	for (int i = 0; i < n_dis-1; i++) {
		if (fehlerAux[i] < twparameter->meshref_kappa*twparameter->meshref_tol &&
		    fehlerAux[i+1] < twparameter->meshref_kappa*twparameter->meshref_tol) {
			//anzI--;
			//I[i]--;
		}

	}
	*/

	// falls durch andere Verfahren bereits Punkte eingefuegt wurden
	for (int i = 0; i < n_dis-1; i++) {
		anzI += I[i];
	}

	if (anzI < M) {
		while (true) {
			if (!FEHLER.empty()){
				fehlerAux = FEHLER.back();
			}

			// Intervall mit groesstem Fehler bestimmen
			max = fehlerAux[0] * pow((1.0 / (1.0 + I[0])), ordVerf-redu+1.0);

			int alpha = 0; // Intervall mit groesstem Fehler

			for (int i = 1; i < n_dis - 1; i++) {
				double aux = fehlerAux[i] * pow((1.0 / (1.0 + I[i])), ordVerf-redu+1.0);
				if (aux > max) {
					max = aux;
					alpha = i;
				}
			}

			// fuer das Abbruchkriterium
			if (I[alpha] == 0 && max < tol) {
				flag = true;
			}
			if (I[alpha] != 0 && max < kappa*tol) {
				flag = true;
			}

			//Abbruchkriterium
			if (anzI > std::min(static_cast<double>(M1), kappa*M)) {
				if (flag || I[alpha] == M1 || anzI == M) {
					break;
				}
			}

			// Vormerkung fuer Punkt
			I[alpha]++;
			anzI++;

			/*
			cout << "........." << endl;
			for (int i = 0; i < n_dis - 1; i++) {
			cout << I[i] << endl;
			}
			cout << "........." << endl;
			*/
		}
	}


	//Ausgabe wie viele Punkte wo eingefuegt werden
	/*
	cout << "I " << anzI << endl;
	for (int i = 0; i < n_dis - 1; i++) {
	cout << I[i] << endl;
	}
	*/
}

vector<double> FullDisTransWorhp::bueskens() {

	vector<double> T_sprungstellen;

	for (int i = 1; i < n_dis; i++) {

		for (int c = 0; c < n_ctrl; c++) {

			/* Multiplikatoren */
			double multi0 = Lambda[u_index(i - 2, c)];
			double multi1 = Lambda[u_index(i - 1, c)];
			double multi2 = Lambda[u_index(i, c)];
			double multi3 = Lambda[u_index(i + 1, c)];

			//cerr << i << ": " << multi1 << ", " << multi2 << ", " << multi3 << endl;

			// Auf und Absprungpunkt
			//if (multi1 > 0 && multi2 > 0 && std::abs(multi3) < 10e-9) {
			if (std::abs(multi2) > 10e-9 && std::abs(multi3) < 10e-9) {
				/* Hilfsvariable */
				double aux = 0;

				if ( i != 1 ) {
					aux = (T[i - 2] * multi1 - T[i - 1] * multi0) / (multi1 - multi0);
					if (aux > T[0]) {
						T_sprungstellen.push_back(aux);
						std::cerr << "++++++++++++++++++++++++++++++++1+ " << aux << endl;
						std::cerr << "T: " << T[i-2] << " " << T[i-1] << " " << T[i] << " " << T[i+1] << endl;
						std::cerr << "M: " << multi0 << " " << multi1 << " " << multi2 << " " << multi3 << endl;
					}
				}
				else {
					aux = (T[i - 1] * multi2 - T[i] * multi1) / (multi2 - multi1);
					if (aux > T[0]) {
						T_sprungstellen.push_back(aux);
						std::cerr << "++++++++++++++++++++++++++++++++2+ " << aux << endl;
						std::cerr << T[i-2] << " " << T[i-1] << " " << T[i] << " " << T[i+1] << endl;
						std::cerr << multi0 << " " << multi1 << " " << multi2 << " " << multi3 << endl;
					}
				}
			}
		}
	}

	return T_sprungstellen;
}


/*
exportiert Lambda in Matlab Datei
TODO: gleiches fuer Mu
*/
void FullDisTransWorhp::ToMATLAB_LambdaMu(const std::string& filename) {

	std::ofstream of(filename);

	of.setf(std::ios::scientific);
	of << std::setprecision(9);

	of << "% TransWORHP-Result for Lambda and Mu" << endl;

	of << "% ";
	for (int i=0; i<n_param; i++) {
		of << std::setw(20) << p(i);
	}
	of << endl;

	//cout << "TERMINATE" << n_dis << endl;
	for (int i=0; i<n_dis; i++) {

		if (phase->freetime) {
			of << std::setw(20) << (T[i] * p(0) );
		} else {
			of << std::setw(20) << T[i];
		}

		for (int j=0; j<n_ode; j++) {
			of << std::setw(20) << Lambda[x_index(i, j)];
		}
		for (int j=0; j<n_ctrl; j++) {
			of << std::setw(20) << Lambda[u_index(i, j)];
		}
		of << endl;
	}
}

/*
liest Lambda aus Matlab Datei und interpoliert Zwischenwerte
TODO: gleiches fuer Mu
*/
void FullDisTransWorhp::FromMATLAB_LambdaMu(const std::string& filename) {

	std::ifstream of(filename);
	string line;

	/// Header
	getline(of, line);
//cout << ":" << line << endl;
	/// Parameter p
	getline(of, line);
	string line2 = string(line,2);
	vector<double> v = ToDoubleArray(line2);
	for (int i=0; i<n_param; i++) {
		X[p_index(i)] = v[i];
	}

	/// Data
	getline(of, line);
//cout << ":" << line << endl;
	vector<double> v_last;

	v = ToDoubleArray(line);
	if (phase->freetime) v[0] = v[0]/p(0);

	v_last = v;
//cout << ":" << v_last << " " << n_dis << endl;
	for (int i=0; i<n_dis; i++) {

		while (T[i] > v[0]) {
//cout << i << "   " << T[i] << " " << v[0] << endl;
			getline(of, line);

			if (!of) {
				v_last = v;
				break;
			} else {
				v_last = v;

				v = ToDoubleArray(line);
				if (phase->freetime) v[0] = v[0]/p(0);
			}
		}

		double vv=1;

		if (v[0]!=v_last[0])
			vv = (T[i] - v_last[0]) / (v[0] - v_last[0]);
//cout << i << "   vv " << vv << endl;

		for (int j=0; j<n_ode; j++) {
			Lambda[x_index(i,j)] = v_last[j+1] + (v[j+1]-v_last[j+1]) * vv;
			/*
			X[x_index(i,j)] = v_last[j+1] + (v[j+1]-v_last[j+1]) * vv;

			if (i>0 && twdiscretization->type==TWdiscretizationType::HermiteSimpson) {
				X[x_index__(2*i-1,j)] = (X[x_index__(2*(i-1),j)] + X[x_index__(2*i,j)])/2;
			}
			*/
		}

		for (int j=0; j<n_ctrl; j++) {

			Lambda[u_index(i,j)] = v_last[j+n_ode+1] + (v[j+n_ode+1]-v_last[j+n_ode+1]) * vv;
			/*
			X[u_index(i,j)] = v_last[j+n_ode+1] + (v[j+n_ode+1]-v_last[j+n_ode+1]) * vv;

			if (i>0 && twdiscretization->type==TWdiscretizationType::HermiteSimpson) {
				X[u_index__(2*i-1,j)] = (X[u_index__(2*(i-1),j)] + X[u_index__(2*i,j)])/2;
			}
			*/
		}


	}

#ifdef TRANSWORHP_GRAPHICS
	if (viewer) viewer->autoScale();
#endif
}

//TODO Mu
exportTW FullDisTransWorhp::outTW() {

	exportTW out;

	// Groesse festlegen
	out.T.reserve(n_dis);
	out.p.reserve(n_param);
	out.x.reserve(twdiscretization->stuetzstellen(n_dis));
	out.u.reserve(twdiscretization->stuetzstellen(n_dis));

	// Diskretisierung
	out.twdiscretisation = twdiscretization;

	// Zeit
	/*if (freetime) {
		for (int i = 0; i < n_dis; i++) {
			out.T.push_back( T[i]*p(0) );
		}
	} else { */
		for (int i = 0; i < n_dis; i++) {
			out.T.push_back( T[i] );
		}
	//}

	// Parameter
	for (int i = 0; i < n_param; i++) {
		out.p.push_back( p(i) );
	}

	// Hilfsvariablen
	vector<double> auxX;
	vector<double> auxU;
	vector<double> auxL;
	//vector<double> auxM;
	auxX.reserve(n_ode);
	auxU.reserve(n_ctrl);
	auxL.reserve(n_var);
	//auxM.reserve(n_con);

	for (int i = 0; i < twdiscretization->stuetzstellen(n_dis); i++) {

		//Zustaende
		for (int j = 0; j < n_ode; j++) {
			auxX.push_back( x__(i, j) );
			auxL.push_back( Lambda[x_index__(i,j)] );
		}
		out.x.push_back(auxX);
		auxX.clear();

		//Steuerung
		for (int j = 0; j < n_ctrl; j++) {
			auxU.push_back( u__(i, j) );
			auxL.push_back( Lambda[u_index__(i,j)] );
		}
		out.u.push_back(auxU);
		auxU.clear();

		out.La.push_back(auxL);
		auxL.clear();
	}

	out.ok = true;

	//cout << "size of output:              " << sizeof(out) << endl;

	return out;
}

// TODO Mu
// TODO Verhalten, wenn Punkte entfernt werden
void FullDisTransWorhp::inTW(exportTW eTW) {

	// Abbruch, falls eTW leer
	if (!eTW.ok) {
		cout << "eTW ist fehlerhaft oder leer" << endl;
		return;
	}

	// Parameter setzen
	for (int i = 0; i < n_param; i++) {
		X[p_index(i)] = eTW.p[i];
	}

	// Zustaende, Steuerungen, Multiplikatoren interpolieren
	if (twdiscretization->type == TWdiscretizationType::HermiteSimpson) { //(quadratisch interpolieren)

		// ersten Eintrag setzen
		for (int k = 0; k < n_ode; k++) {
			X[x_index(0,k)] = eTW.x[0][k];
		}
		for (int k = 0; k < n_ctrl; k++) {
			X[u_index(0,k)] = eTW.u[0][k];
		}

		for (int i = 1, j = 1; i < n_dis; j++) {

			while (T[i] <= eTW.T[j] && i < n_dis) {

				double x0 = eTW.T[j-1];
				double h = eTW.T[j]-x0;

				// Zustand
				for (int k = 0; k < n_ode; k++) {

					double f0 = eTW.x[2*j-2][k];
					double f1 = eTW.x[2*j-1][k];
					double f2 = eTW.x[2*j][k];

					// an Zwischenstelle
					double TT = (T[i]+T[i-1])/2.;
					X[x_index__(2*i-1,k)] = 2.0*f0*(TT-h/2.0-x0)/(h*h)*(TT-h-x0)-4.0*f1*(TT-x0)/(h*h)*(TT-h-x0)+2.0*f2*(TT-x0)/(h*h)*(TT-h/2.0-x0);

					// an Gitterpunkt
					X[x_index(i,k)] = 2.0*f0*(T[i]-h/2.0-x0)/(h*h)*(T[i]-h-x0)-4.0*f1*(T[i]-x0)/(h*h)*(T[i]-h-x0)+2.0*f2*(T[i]-x0)/(h*h)*(T[i]-h/2.0-x0);

					// Lambda ///////////

					f0 = eTW.La[2*j-2][k];
					f1 = eTW.La[2*j-1][k];
					f2 = eTW.La[2*j][k];

					// an Zwischenstelle
					TT = (T[i]+T[i-1])/2.;
					Lambda[x_index__(2*i-1,k)] = 2.0*f0*(TT-h/2.0-x0)/(h*h)*(TT-h-x0)-4.0*f1*(TT-x0)/(h*h)*(TT-h-x0)+2.0*f2*(TT-x0)/(h*h)*(TT-h/2.0-x0);

					// an Gitterpunkt
					Lambda[x_index__(2*i-1,k)] = 2.0*f0*(T[i]-h/2.0-x0)/(h*h)*(T[i]-h-x0)-4.0*f1*(T[i]-x0)/(h*h)*(T[i]-h-x0)+2.0*f2*(T[i]-x0)/(h*h)*(T[i]-h/2.0-x0);
				}

				// Steuerung
				for (int k = 0; k < n_ctrl; k++) {

					double f0 = eTW.u[2*j-2][k];
					double f1 = eTW.u[2*j-1][k];
					double f2 = eTW.u[2*j][k];

					// an Zwischenstelle
					double TT = (T[i]+T[i-1])/2.;
					X[u_index__(2*i-1,k)] = 2.0*f0*(TT-h/2.0-x0)/(h*h)*(TT-h-x0)-4.0*f1*(TT-x0)/(h*h)*(TT-h-x0)+2.0*f2*(TT-x0)/(h*h)*(TT-h/2.0-x0);

					// an Gitterpunkt
					X[u_index(i,k)] = 2.0*f0*(T[i]-h/2.0-x0)/(h*h)*(T[i]-h-x0)-4.0*f1*(T[i]-x0)/(h*h)*(T[i]-h-x0)+2.0*f2*(T[i]-x0)/(h*h)*(T[i]-h/2.0-x0);

					// Lambda ///////////

					f0 = eTW.La[2*j-2][k+n_ode];
					f1 = eTW.La[2*j-1][k+n_ode];
					f2 = eTW.La[2*j][k+n_ode];

					// an Zwischenstelle
					TT = (T[i]+T[i-1])/2.;
					Lambda[u_index__(2*i-1,k)] = 2.0*f0*(TT-h/2.0-x0)/(h*h)*(TT-h-x0)-4.0*f1*(TT-x0)/(h*h)*(TT-h-x0)+2.0*f2*(TT-x0)/(h*h)*(TT-h/2.0-x0);

					// an Gitterpunkt
					Lambda[u_index__(2*i-1,k)] = 2.0*f0*(T[i]-h/2.0-x0)/(h*h)*(T[i]-h-x0)-4.0*f1*(T[i]-x0)/(h*h)*(T[i]-h-x0)+2.0*f2*(T[i]-x0)/(h*h)*(T[i]-h/2.0-x0);
				}

				i++;
			}
		}

	} else { // Euler und Trapez (linear interpolieren)

		for (int i = 0, j = 0; i < n_dis; i++, j++) {

			while (T[i] < eTW.T[j]) {

				for (int k = 0; k < n_ode; k++) {
					X[x_index(i,k)] = eTW.x[j-1][k] + (eTW.x[j][k]-eTW.x[j-1][k])/(eTW.T[j]-eTW.T[j-1])*(T[i]-eTW.T[j-1]);

					Lambda[x_index(i,k)] = eTW.La[j-1][k] + (eTW.La[j][k]-eTW.La[j-1][k])/(eTW.T[j]-eTW.T[j-1])*(T[i]-eTW.T[j-1]);
				}
				for (int k = 0; k < n_ctrl; k++) {
					X[u_index(i,k)] = eTW.u[j-1][k] + (eTW.u[j][k]-eTW.u[j-1][k])/(eTW.T[j]-eTW.T[j-1])*(T[i]-eTW.T[j-1]);

					Lambda[u_index(i,k)] = eTW.La[j-1][k+n_ode] + (eTW.La[j][k+n_ode]-eTW.La[j-1][k+n_ode])/(eTW.T[j]-eTW.T[j-1])*(T[i]-eTW.T[j-1]);
				}

				i++;
			}

			if (T[i] == eTW.T[j]) {
				for (int k = 0; k < n_ode; k++) {
					X[x_index(i,k)] = eTW.x[j][k];
					Lambda[x_index(i,k)] = eTW.La[j][k];
				}
				for (int k = 0; k < n_ctrl; k++) {
					X[u_index(i,k)] = eTW.u[j][k];
					Lambda[u_index(i,k)] = eTW.La[j][k+n_ode];
				}
			}
		}
	}

#ifdef TRANSWORHP_GRAPHICS
	if (viewer) viewer->autoScale();
#endif

}

}
