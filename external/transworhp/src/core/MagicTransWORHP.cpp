#ifdef WIN32
#include <windows.h>
#endif

#include "MagicTransWORHP.h"
#include <iomanip>
#include "worhp_info.h"
#include "conversion.h"
#include "xmlio.h"
#include <sstream>
#include "../base/twstatus.h"

using namespace std;

MagicTransWorhp::MagicTransWorhp(const char *s, int dis, int ode, int ctrl, int param, int rand, int neben)
	: TransWorhp(s,dis,ode,ctrl,param,rand,neben), Lagrangian(0) {

	MagicDouble::Derivate = 1;

}


void MagicTransWorhp::localinit() {
	magic_X1 = MagicDouble::NEW(n_ode);
	magic_U1 = MagicDouble::NEW(n_ctrl);
	magic_P  = MagicDouble::NEW(n_param);

	nomagic_DX1 = MagicDouble::NO_DEP_NEW(n_ode);
	nomagic_DX12= MagicDouble::NO_DEP_NEW(n_ode);
	nomagic_DX2 = MagicDouble::NO_DEP_NEW(n_ode);

	nomagic_X1 = MagicDouble::NO_DEP_NEW(n_ode);
	nomagic_X12= MagicDouble::NO_DEP_NEW(n_ode);
	nomagic_X2 = MagicDouble::NO_DEP_NEW(n_ode);

	nomagic_U1 = MagicDouble::NO_DEP_NEW(n_ctrl);
	nomagic_U12= MagicDouble::NO_DEP_NEW(n_ctrl);
	nomagic_U2 = MagicDouble::NO_DEP_NEW(n_ctrl);

	nomagic_P = MagicDouble::NO_DEP_NEW(n_param);

	magic_XX = MagicDouble::NEW( n_var );
	nomagic_XX = MagicDouble::NO_DEP_NEW( n_var );
	nomagic_GG = MagicDouble::NO_DEP_NEW( n_con );

}

MagicTransWorhp::~MagicTransWorhp() {

	delete []magic_X1;
	delete []magic_U1;
	delete []magic_P;

	delete []nomagic_DX1;
	delete []nomagic_DX12;
	delete []nomagic_DX2;

	delete []nomagic_X1;
	delete []nomagic_X12;
	delete []nomagic_X2;

	delete []nomagic_U1;
	delete []nomagic_U12;
	delete []nomagic_U2;

	delete []nomagic_P;

	delete []magic_XX;
	delete []nomagic_XX;
	delete []nomagic_GG;

}


double MagicTransWorhp::obj() {

	MagicDouble::COPY(nomagic_XX,X,n_var);
 	MagicDouble F(0);

	obj(&F, nomagic_XX);

	return F.F();

}

bool MagicTransWorhp::obj_structure(DiffStructure &s) {

 	MagicDouble F(0);

	obj(&F, magic_XX);
	for (int k=0; k<n_var; k++) {
		const double *z = F.DF(&magic_XX[k]);
		if (z)
			s(0, k);
	}

	return true;
}


void MagicTransWorhp::DF_calculate(WorhpMatrix &DF, double ScaleObj) {

	MagicDouble::COPY(magic_XX,X,n_var);
	MagicDouble F(0);

	obj(&F, magic_XX);

	for (int k=0; k<DF.nnz; k++) {
		const double *z = F.DF(&magic_XX[DF.row[k]-1]);
		if (z)
		DF.val[k] = ScaleObj * *z;
	}

	twcount_calls.obj.df++;

}


void MagicTransWorhp::ode(double *dx, double t, const double *x, const double *u, const double *p) {

	MagicDouble::COPY(nomagic_X1,x,n_ode);
	MagicDouble::COPY(nomagic_U1,u,n_ctrl);
	MagicDouble::COPY(nomagic_P,p,n_param);
	MagicDouble::COPY(nomagic_DX1,dx,n_ode);

	ode(nomagic_DX1,t,nomagic_X1,nomagic_U1,nomagic_P);

	for (int i=0;i<n_ode;i++) {
		dx[i] = nomagic_DX1[i].F();
	}

}

bool MagicTransWorhp::ode_structure(DiffStructure &/*s*/) {

	double t=1;

	ode(nomagic_DX1,t, //T[dis],
	    magic_X1,  // Zustand X
	    magic_U1, // Steuerung U
	    magic_P);  // Parameter P

	for (int j=0; j<n_ode; j++) {

		for (int k=0;k<n_ode;k++) {
			const double *z = nomagic_DX1[j].DF(&magic_X1[k]);
			if (z) DS_ode(j, x_indexode(k));
		}

		for (int k=0;k<n_ctrl;k++)  {
			const double *z = nomagic_DX1[j].DF(&magic_U1[k]);
			if (z) DS_ode(j, u_indexode(k));
		}

		for (int k=0;k<n_param;k++)  {
			const double *z = nomagic_DX1[j].DF(&magic_P[k]);
			if (z) DS_ode(j, p_indexode(k));
		}

	}

	return true;
}



void MagicTransWorhp::RechteSeite(double */*G1*/, double /*t*/, int /*dis*/) {

	cout << "ERROR! RechteSeite() macht keinen Sinn!!!!" << endl;

}



void MagicTransWorhp::Trapez(double* G1, int dis) {

	MagicDouble::COPY(nomagic_X1,&X[(n_ode+n_ctrl)* twdiscretization->punkte()*dis],n_ode);
	MagicDouble::COPY(nomagic_X2,&X[(n_ode+n_ctrl)* twdiscretization->punkte()*(dis+1)],n_ode);

	MagicDouble::COPY(nomagic_U1,&X[(n_ode+n_ctrl)* twdiscretization->punkte()*dis + n_ode],n_ctrl);
	MagicDouble::COPY(nomagic_U2,&X[(n_ode+n_ctrl)* twdiscretization->punkte()*(dis+1) + n_ode],n_ctrl);

	MagicDouble::COPY(nomagic_P,&X[(n_ode+n_ctrl)* twdiscretization->stuetzstellen(n_dis)],n_param);

	double T1=T[dis];
	double T2=T[dis+1];

	twcount_calls.ode.f +=2;
	ode(nomagic_DX1,T1,
	    nomagic_X1,  // Zustand X
	    nomagic_U1, // Steuerung U
	    nomagic_P);  // Parameter P
	ode(nomagic_DX2,T2,
	    nomagic_X2,  // Zustand X
	    nomagic_U2, // Steuerung U
	    nomagic_P);  // Parameter P

	for (int j=0; j<n_ode; j++) {
		G1[j] = x(dis+1,j) - x(dis,j) - (T[dis+1] - T[dis]) * (nomagic_DX1[j].F()+nomagic_DX2[j].F())/2;
	}

}


void MagicTransWorhp::HermiteSimpson(double* G1, double *G2, int dis) {

	MagicDouble::COPY(nomagic_X1,&X[(n_ode+n_ctrl)* twdiscretization->punkte()*dis],n_ode);
	MagicDouble::COPY(nomagic_X12,&X[(n_ode+n_ctrl)* (twdiscretization->punkte()*dis+1)],n_ode);
	MagicDouble::COPY(nomagic_X2,&X[(n_ode+n_ctrl)* twdiscretization->punkte()*(dis+1)],n_ode);

	MagicDouble::COPY(nomagic_U1,&X[(n_ode+n_ctrl)* twdiscretization->punkte()*dis + n_ode],n_ctrl);
	MagicDouble::COPY(nomagic_U12,&X[(n_ode+n_ctrl)* (twdiscretization->punkte()*dis+1) + n_ode],n_ctrl);
	MagicDouble::COPY(nomagic_U2,&X[(n_ode+n_ctrl)* twdiscretization->punkte()*(dis+1) + n_ode],n_ctrl);

	MagicDouble::COPY(nomagic_P,&X[(n_ode+n_ctrl)* twdiscretization->stuetzstellen(n_dis)],n_param);

	double T1=T[dis];
	double T2=T[dis+1];
	double T12=(T1+T2)*.5;

	twcount_calls.ode.f +=3;
	ode(nomagic_DX1,T1,
	    nomagic_X1,  // Zustand X
	    nomagic_U1, // Steuerung U
	    nomagic_P);  // Parameter P
	ode(nomagic_DX2,T2,
	    nomagic_X2,  // Zustand X
	    nomagic_U2, // Steuerung U
	    nomagic_P);  // Parameter P
	ode(nomagic_DX12,T12,
	    nomagic_X12,  // Zustand X
	    nomagic_U12, // Steuerung U
	    nomagic_P);  // Parameter P

	// Hermite
	for (int j=0; j<n_ode; j++) {

		//MagicDouble *x1 = nomagic_X1 + j;
		//MagicDouble *x2 = nomagic_X2 + j;
		//MagicDouble *x12 = nomagic_X12 + j;

		G1[j] = X[(n_ode + n_ctrl) * (2 * dis+1) + j]
		        - .5*(x(dis+1,j) + x(dis,j))
		        - .125*(T[dis+1] - T[dis]) * (nomagic_DX1[j].F()-nomagic_DX2[j].F());

	}
	// Simpson
	for (int j=0; j<n_ode; j++) {
		//MagicDouble *x1 = nomagic_X1 + j;
		//MagicDouble *x2 = nomagic_X2 + j;

		G2[j] = x(dis+1,j) - x(dis,j) - (T[dis+1] - T[dis])/6. * (nomagic_DX1[j].F()+nomagic_DX2[j].F()+4.*nomagic_DX12[j].F());

	}

}




// Trapez mit 1. + 2. Ableitung
void MagicTransWorhp::Trapez(MagicDouble* XX, MagicDouble* G1, int dis) {

	MagicDouble *X1 = &XX[(n_ode+n_ctrl)* twdiscretization->punkte()*dis];
	MagicDouble *X2 = &XX[(n_ode+n_ctrl)* twdiscretization->punkte()*(dis+1)];

	MagicDouble *U1 = &XX[(n_ode+n_ctrl)* twdiscretization->punkte()*dis + n_ode];
	MagicDouble *U2 = &XX[(n_ode+n_ctrl)* twdiscretization->punkte()*(dis+1) + n_ode];

	MagicDouble *P  = &XX[(n_ode+n_ctrl)* twdiscretization->stuetzstellen(n_dis)];

	double T1=T[dis];
	double T2=T[dis+1];

	twcount_calls.ode.f +=2;
	ode(nomagic_DX1,T1,
	    X1,  // Zustand X
	    U1, // Steuerung U
	    P);  // Parameter P
	ode(nomagic_DX2,T2,
	    X2,  // Zustand X
	    U2, // Steuerung U
	    P);  // Parameter P


	for (int j=0; j<n_ode; j++) {
		G1[j] = XX[x_index(dis+1,j)] - XX[x_index(dis,j)] - ((T2 - T1)/2)*(nomagic_DX1[j]+nomagic_DX2[j]);
		//G1[j] = x(dis+1,j) - x(dis,j) - (T[dis+1] - T[dis]) * (nomagic_DX1[j]+nomagic_DX2[j])/2;
	}


}


// Hermite-Simpson mit 1. + 2. Ableitung
void MagicTransWorhp::HermiteSimpson(MagicDouble* XX, MagicDouble* G1, MagicDouble *G2, int dis) {

	MagicDouble *X1 = &XX[(n_ode+n_ctrl)* twdiscretization->punkte()*dis];
	MagicDouble *X12 = &XX[(n_ode+n_ctrl)* (twdiscretization->punkte()*dis+1)];
	MagicDouble *X2 = &XX[(n_ode+n_ctrl)* twdiscretization->punkte()*(dis+1)];

	MagicDouble *U1 = &XX[(n_ode+n_ctrl)* twdiscretization->punkte()*dis + n_ode];
	MagicDouble *U12 = &XX[(n_ode+n_ctrl)* (twdiscretization->punkte()*dis+1) + n_ode];
	MagicDouble *U2 = &XX[(n_ode+n_ctrl)* twdiscretization->punkte()*(dis+1) + n_ode];

	MagicDouble *P  = &XX[(n_ode+n_ctrl)* twdiscretization->stuetzstellen(n_dis)];

	double T1=T[dis];
	double T2=T[dis+1];
	double T12=(T[dis]+T[dis+1])*.5;

	twcount_calls.ode.f +=3;
	ode(nomagic_DX1,T1,
	    X1,  // Zustand X
	    U1, // Steuerung U
	    P);  // Parameter P
	ode(nomagic_DX2,T2,
	    X2,  // Zustand X
	    U2, // Steuerung U
	    P);  // Parameter P
	ode(nomagic_DX12,T12,
	    X12,  // Zustand X
	    U12, // Steuerung U
	    P);  // Parameter P

	// Hermite
	for (int j=0; j<n_ode; j++) {

		MagicDouble *x1 = X1 + j;
		MagicDouble *x2 = X2 + j;
		MagicDouble *x12 = X12 + j;
		G1[j] = *x12
		        - .5*(*x1+*x2)
		        -( .125*(T[dis+1] - T[dis])) * (nomagic_DX1[j]-nomagic_DX2[j]);
	}



	// Simpson
	for (int j=0; j<n_ode; j++) {

		MagicDouble *x1 = X1 + j;
		MagicDouble *x2 = X2 + j;

		G2[j] = *x2 - *x1 - ((T[dis+1] - T[dis])/6.) * (nomagic_DX1[j]+nomagic_DX2[j]+4.*nomagic_DX12[j]);
	}


/*	for (int j=0; j<n_ode; j++) {
		G1[j] = X[(n_ode + n_ctrl) * (2 * dis+1) + j]
		        - .5*(x(dis+1,j) + x(dis,j))
		        - .125*(T[dis+1] - T[dis]) * (nomagic_DX1[j].F()-nomagic_DX2[j].F());
	}

	// Simpson
	for (int j=0; j<n_ode; j++) {
		G2[j] = x(dis+1,j) - x(dis,j) - (T[dis+1] - T[dis])/6. * (nomagic_DX1[j].F()+nomagic_DX2[j].F()+4.*nomagic_DX12[j].F());
	}*/
}



void MagicTransWorhp::DG_diff_ode(double t, int dis, int active_index) {

	//TransWorhp::DG_diff_ode(t,dis,active_index);
	//return;


	// in DS_ode ablegen
	// dis von 0 bis Disc.stuetzstellen(n_dis)

	DS_ode.activeindex = active_index;

	// Analytisch!!!
	MagicDouble::COPY(magic_X1,&X[(n_ode+n_ctrl)*dis],n_ode);
	MagicDouble::COPY(magic_U1,&X[(n_ode+n_ctrl)*dis + n_ode],n_ctrl);
	MagicDouble::COPY(nomagic_P,&X[(n_ode+n_ctrl)* twdiscretization->stuetzstellen(n_dis)],n_param);

	ode(nomagic_DX1,t, //T[dis],
	    magic_X1,  // Zustand X
	    magic_U1, // Steuerung U
	    nomagic_P);  // Parameter P


	for (int j=0;j<DS_ode.n_eq;j++) {
		for (int k=0;k<DS_ode.n_diff;k++) {

			if (DS_ode.check(j,k)) {

				if (k<n_ode) {
				const double *z = nomagic_DX1[j].DF(&magic_X1[k]);
				if (z) DS_ode(j, k) = *z;
				//else cout << "WERWER" << j << " " << k << endl;
				}
				else if (k<n_ode+n_ctrl) {
				const double *z = nomagic_DX1[j].DF(&magic_U1[k-n_ode]);
				if (z) DS_ode(j, k) = *z;
				//else cout << "WERWER u  " << j << " " << k << endl;
				}
				/*else  {
				const double *z = nomagic_DX1[j].DF(&magic_P[k-n_ode-n_ctrl]);
				if (z) DS_ode(j, k) = *z;
				else cout << "WERWER p " << j << " " << k<< endl;
				}*/
			}
		}
	}

	/*for (int j=0; j<n_ode; j++) {

		for (int k=0;k<n_ode;k++) {
			const double *z = nomagic_DX1[j].DF(&magic_X1[k]);
			if (z) DS_ode(j, x_indexode(k)) = *z;
		}

		for (int k=0;k<n_ctrl;k++)  {
			const double *z = nomagic_DX1[j].DF(&magic_U1[k]);
			if (z) DS_ode(j, u_indexode(k)) = *z;
		}

		for (int k=0;k<n_param;k++)  {
			const double *z = nomagic_DX1[j].DF(&nomagic_P[k]);
			if (z) DS_ode(j, p_indexode(k)) = *z;
		}

	}*/

	twcount_calls.ode.df++;
}


void MagicTransWorhp::DG_diff_ode_p(double t, int l, int dis, int active_index) {

	// in DS_ode ablegen
	// dis von 0 bis Disc.stuetzstellen(n_dis)

	DS_ode.activeindex = active_index;

	// Analytisch!!!
	MagicDouble::COPY(nomagic_X1,&X[(n_ode+n_ctrl)* dis],n_ode);
	MagicDouble::COPY(nomagic_U1,&X[(n_ode+n_ctrl)* dis + n_ode],n_ctrl);
	MagicDouble::COPY(magic_P,&X[(n_ode+n_ctrl)* twdiscretization->stuetzstellen(n_dis)],n_param);



	ode(nomagic_DX1,t, //T[dis],
	    nomagic_X1,  // Zustand X
	    nomagic_U1, // Steuerung U
	    magic_P);  // Parameter P

	/*for (int j=0;j<n_ode;j++)  {
		const double *z = nomagic_DX1[j].DF(&magic_P[l]);
		if (z) DS_ode(j, p_indexode(l)) = *z;
	}*/

	int k = p_indexode(l);

	for (int j=0;j<DS_ode.n_eq;j++) {
		//for (int k=0;k<DS_ode.n_diff;k++) {

			if (DS_ode.check(j,k)) {

				/*if (k<n_ode) {
				const double *z = nomagic_DX1[j].DF(&magic_X1[k]);
				if (z) DS_ode(j, k) = *z;
				else cout << "WERWER" << j << " " << k << endl;
				}
				else if (k<n_ode+n_ctrl) {
				const double *z = nomagic_DX1[j].DF(&magic_U1[k-n_ode]);
				if (z) DS_ode(j, k) = *z;
				else cout << "WERWER u  " << j << " " << k << endl;
				}*/

				const double *z = nomagic_DX1[j].DF(&magic_P[k-n_ode-n_ctrl]);
				if (z) DS_ode(j, k) = *z;
				//else cout << "WERWER p " << j << " " << k<< endl;

			}
		}



	twcount_calls.ode.dfp++;

}


void MagicTransWorhp::infoMagic() {

	MyStatus("TransWORHP", "Using Automatic Differentiation for ode system", Status::WARN);

}






void MagicTransWorhp::Objective(MagicDouble *F, double ScaleObj) {

	//MagicDouble::DEBUG(magic_XX,n_var);

	//MagicDouble::COPY(nomagic_XX,X,n_var);
	obj(F, magic_XX);

	*F = *F * ScaleObj;

}




void MagicTransWorhp::Constraints3(MagicDouble* XX, MagicDouble *GG) {

	if (twdiscretization->type==TWdiscretizationType::Trapez) {

		for (int dis=0; dis<n_dis-1; dis++) {
			//cout << "Constraints3    dis" << dis << endl;
			MagicDouble *G1 = &GG[dis*n_ode];
			Trapez(XX, G1, dis);
		}
	} else if (twdiscretization->type==TWdiscretizationType::HermiteSimpson) {

		for (int dis=0; dis<n_dis-1; dis++) {

			MagicDouble *G1 = &GG[2*dis*n_ode];
			MagicDouble *G2 = &GG[(2*dis+1)*n_ode];
			HermiteSimpson(XX, G1, G2, dis);


		}
	}

/*	if (n_rand) {
		count.rand.f++;
		rand( &GG[Disc.stufen * (n_dis-1) * n_ode] );
	}

	if (n_neben) {
		for (int dis=0; dis<n_dis; dis++) {
			double *G1 = &GG[Disc.stufen * (n_dis-1) * n_ode + n_rand + dis*n_neben];

			double *X1 = &X[(n_ode+n_ctrl)* Disc.punkte() * dis];
			double *U1 = &X[(n_ode+n_ctrl)* Disc.punkte() * dis + n_ode];
			double *P = &X[(n_ode+n_ctrl)* Disc.stuetzstellen(n_dis)];
			count.neben.f++;

			neben(G1,T[dis],X1,U1,P);
		}
	}*/

}

#include "../base/progressbar.hpp"


int MagicTransWorhp::HM_structure(int hessianstructure, WorhpMatrix &/*DF*/, WorhpMatrix &/*DG*/, WorhpMatrix *HM) {


	MagicDouble::Derivate=2;

	int ind = 0;


	if (hessianstructure==4) { // To slow --- 5 verwenden!


		ProgressBar P(n_var*(n_var-1)/2);
		int pp=0;


		MagicDouble F(0);
		//MagicDouble::COPY(nomagic_XX,worhp_o->X,n_var);

		for (int y = 0; y < n_var; ++y) {
			for (int x = y+1; x < n_var; ++x) {

				pp++;
				P.status(pp);

				//cout << y << " " << x << endl;

				nomagic_XX[x].link();
				nomagic_XX[y].link();

				//MagicDouble::DEBUG(nomagic_XX,n_var);

				//	cout << "Eval" << endl;

				obj(&F, nomagic_XX);
				//	cout << "Eval0" << endl;
				Constraints3(nomagic_XX, nomagic_GG);


				//MagicDouble::DEBUG(nomagic_GG,n_con);
				//	cout << "ok" << endl;

				nomagic_XX[x].unlink();
				nomagic_XX[y].unlink();


			//	double lxx = 0;

			//	const double *z = F.DDF(&nomagic_XX[x-1], &nomagic_XX[y-1]);

			//	if (z) lxx=*z;

				for (int i=0;i<n_con;i++) {
					F = F  + /*Mu[i] **/  nomagic_GG[i];
				}

				const double *z = F.DDF(&nomagic_XX[x], &nomagic_XX[y]);
				if (z) {
					//cout << "ZZZ" << endl;
					if (HM) HM->row[ind] = x + 1;
					if (HM) HM->col[ind] = y + 1;
					ind++;
				}

			}
		}

	}

	else if (hessianstructure==5) {

		cout << "GO HM" << endl;
		if (HM==nullptr) {
			obj(&Lagrangian, magic_XX);
			Constraints3(magic_XX,nomagic_GG);

			for (int i=0;i<n_con;i++) {
				Lagrangian = Lagrangian  + /*Mu[i] **/  nomagic_GG[i];
			}
		}
		cout << "GO HM 2" << endl;

		for (int y = 0; y < n_var; ++y) {
			for (int x = y+1; x < n_var; ++x) {

				const double *z = Lagrangian.DDF(&magic_XX[x], &magic_XX[y]);
				if (z) {
					if (HM) HM->row[ind] = x + 1;
					if (HM) HM->col[ind] = y + 1;
					ind++;
				}
			}
		}
		cout << "GO HM 3" << endl;
	}

	else {
	//	return TransWorhp::HM_structure(DF,DG,HM);
	}


	// Diagonale
	for (int i = 0; i < n_var; ++i) {

		if (HM) HM->row[ind] = i + 1;
		if (HM) HM->col[ind] = i + 1;
		ind++;

	}

MagicDouble::Derivate=1;
	return ind;
}

// wenig Variablen -- oft aufrufen (und unten noch umgekehrt!!!)
void MagicTransWorhp::HM_calculate4(WorhpMatrix &/*DF*/, WorhpMatrix &/*DG*/, WorhpMatrix &HM, double ScaleObj, double *Mu) {

	MagicDouble::Derivate=2;
	//cout << "HM_calculate4" << endl;

	MagicDouble F(0);
	MagicDouble::COPY(nomagic_XX,X,n_var);

	for (int ind=0;ind<HM.nnz;ind++) {

		int x = HM.col[ind];
		int y = HM.row[ind];

		if (x==y) {

			nomagic_XX[x-1].link();
		//	MagicDouble::DEBUG(nomagic_XX,n_var);

			Objective(&F, ScaleObj);
			Constraints3(nomagic_XX, nomagic_GG);

			nomagic_XX[x-1].unlink();

		//	cout << "F: " << F << endl;
		//	MagicDouble::DEBUG(nomagic_GG,n_con);
		//	MagicDouble::DEBUG(nomagic_XX,n_var);
		//	cout << endl;

		}

		else {

			nomagic_XX[x-1].link();
			nomagic_XX[y-1].link();

			Objective(&F, ScaleObj);
			Constraints3(nomagic_XX, nomagic_GG);

			nomagic_XX[x-1].unlink();
			nomagic_XX[y-1].unlink();

			//ddf = (f1-f2-f3+f4) / (4*eps*eps);


		}


		double lxx = 0;

		const double *z = F.DDF(&nomagic_XX[x-1], &nomagic_XX[y-1]);
		if (z) lxx=*z;

		for (int i=0;i<n_con;i++) {
			const double *z = nomagic_GG[i].DDF(&nomagic_XX[x-1], &nomagic_XX[y-1]);
			if (z)
				lxx += Mu[i] * *z;
		}


		HM.val[ind] = lxx;

	}

	MagicDouble::Derivate=1;

}


// viele Variablen -- einmal aufrufen
void MagicTransWorhp::HM_calculate5(WorhpMatrix &/*DF*/, WorhpMatrix &/*DG*/, WorhpMatrix &HM, double ScaleObj, double *Mu) {

	//cout << "HM_calculate5" << endl;
	MagicDouble::Derivate=2;

	MagicDouble F(0);
	MagicDouble::COPY(magic_XX,X,n_var);

	obj(&F, magic_XX);
	F *= ScaleObj;

	Constraints3(magic_XX, nomagic_GG);

	for (int i=0;i<n_con;i++) {
		F += Mu[i] * nomagic_GG[i];
	}
	//cout << F << endl;

	for (int ind=0;ind<HM.nnz;ind++) {

		int x = HM.col[ind];
		int y = HM.row[ind];

		const double *z = F.DDF(&magic_XX[x-1], &magic_XX[y-1]);
		if (z) HM.val[ind] = *z;

		else HM.val[ind] = 0;
	}

	MagicDouble::Derivate=1;

}


void MagicTransWorhp::ParseXML(XMLNode *xmlmain) {

	if (!xmlmain) return;

	TransWorhp::ParseXML(xmlmain);

	XMLNode *n = xmlmain->GetFirstChild("MAGIC_HESSIAN");
	if (n) {
/*
		XMLNode *nn = n->GetFirstChild("STRUCTURE");

		if (nn) {
			string s = nn->GetText();

			hessianstructure = ToInt(s);

			if (s=="Diagonal") {
				hessianstructure = 0;
			}
			if (s=="Full") {
				hessianstructure = 1;
			}
			if (s=="Odeblocks2") {
				hessianstructure = 2;
			}
			if (s=="Odeblocks") {
				hessianstructure = 3;
			}
		}

		nn = n->GetFirstChild("VALUES");

		if (nn) {
			string s = nn->GetText();

	//		hessianvalues = ToInt(s);

			if (s=="DiffDiff") {
	//			hessianvalues = 0;
			}
			if (s=="DiffDG") {
	//			hessianvalues = 1;
			}

		}*/
	}
}

