#include "GaussPmTransWORHP.h"

#include "legendreGaussPoints.h"
#include "legendreGaussWeights.h"

#include "butcher.h"
#include "conversion.h"
#include "TWproblem.h"
#include "TWfolder.h"
#include "twstatus.h"

#ifdef TRANSWORHP_GRAPHICS
#include "Viewer.h"
#endif

#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::stringstream;

namespace tw {

GaussPmTransWorhp::GaussPmTransWorhp(TransWorhpProblem* ph, TWparameter *twparam)
	: PmTransWorhp(ph,twparam)
{

	transworhp_type = TransWORHP_type::pseudospectral_gauss;
	
	
	time.resize(n_dis);
	legendreGaussTime();
	
	/*
	for (auto ele : time) {
		cout << ele << " ";
	}
	exit(1);
	//*/
	
	n_collo = n_dis;
	n_dis = n_collo+2; // Anfangs- und Endwert
	
	T.resize(n_dis);
	for (int i = 0; i < n_dis; i++) {
		T[i] = timeToTzeroTfinal(time[i]);
	}
	
	// Groesse Integral
	lagrange_integral.resize(n_integral*(n_collo+1));
	
	gauss_weights = getLegendreGaussWeights(n_collo);
	
	if (!smoothMode) {
		displayPoints = n_dis;
	}

	displayTime.reserve(displayPoints);
	for (int i = 0; i < displayPoints; i++) {
		const double t = i/(displayPoints-1.)*(phase->tf-phase->t0)+phase->t0;
		//cout << t << endl;
		displayTime.push_back(t);
	}

	displayXU.resize((n_ode+n_ctrl)*displayPoints);

	//Pseudospektral-Matrix fuellen
	fillGaussPseudospectralMatrixBarycentric();
	//fillGaussPseudospectralMatrixDirect();
	
	/* //Pseudospektral-Differential-Matrix
	for (size_t k = 0; k < PDM_D.size(); k++) {
		for (size_t i = 0; i < PDM_D[k].size(); i++) {
			cout << std::setprecision(3) << std::setw(5) << PDM_D[k][i] << "\t";
		}
		cout << endl;
	}
	*/
	
	/*
	for (auto ele : gauss_weights) {
		cout << ele << endl;
	}
	*/
	
	//exit(1);
}


GaussPmTransWorhp::~GaussPmTransWorhp() {

}


std::string GaussPmTransWorhp::type_G(int row) const {

	stringstream a;

	if (row < n_collo*n_ode) {
		a << " COL  " << std::setw(2) << row%n_ode << " ";
	}
	else if (row < n_collo*n_ode+n_ode) {
		a << " COLe " << std::setw(2) << (row-n_collo*n_ode)%n_ode << " ";
	}
	else if (row < n_collo*n_ode+n_ode + n_rand) {
		a << " RAND " << std::setw(2) << (row-n_collo*n_ode+n_ode)%n_rand << " ";
	}
	else if (row < n_collo*n_ode+n_ode+n_rand + n_neben*n_dis) {
		a << " NB   " << std::setw(2) << (row-n_collo*n_ode+n_ode-n_rand-n_neben)%n_neben << " ";
	}

	return a.str();
}


void GaussPmTransWorhp::fromMATLAB_impl(std::ifstream &/*stream*/) {
/* nicht getestet
	std::ifstream of(filename);
	string line;

	/// Header
	std::getline(of, line);
	/// Parameter p
	std::getline(of, line);
	string line2 = string(line,2);
	vector<double> v = ToDoubleArray(line2);
	for (int i=0; i<n_param; i++) {
		X[p_index(i)] = v[i];
	}

	/// Data
	std::getline(of, line);
	vector<double> v_last;

	v = ToDoubleArray(line);
	if (phase->freetime) v[0] = v[0]/p(0);

	v_last = v;
	for (int i = 0; i < n_dis; i++) {

		while (T[i] > v[0]) {
			std::getline(of, line);

			if (!of) {
				v_last = v;
				break;
			} else {
				v_last = v;

				v = ToDoubleArray(line);
				if (phase->freetime) v[0] = v[0]/p(0);
			}
		}

		double vv = 1;

		if (v[0] != v_last[0]) {
			vv = (T[i] - v_last[0]) / (v[0] - v_last[0]);
		}

		for (int j=0; j<n_ode; j++) {
			X[x_index(i,j)] = v_last[j+1] + (v[j+1]-v_last[j+1]) * vv;
		}

		for (int j=0; j<n_ctrl; j++) {
			X[u_index(i,j)] = v_last[j+n_ode+1] + (v[j+n_ode+1]-v_last[j+n_ode+1]) * vv;
		}

	}
*/
#ifdef TRANSWORHP_GRAPHICS
	if (viewer) viewer->autoScale();
#endif
}


void GaussPmTransWorhp::legendreGaussTime() {

	time = getLegendreGaussPoints(n_dis);
	time.push_back(-1.0); // Anfangswert (KEIN Kollokationspunkt!)
	time.push_back(+1.0); // Endwert (KEIN Kollokationspunkt!)

	std::sort(time.begin(),time.end());
}


void GaussPmTransWorhp::fillGaussPseudospectralMatrixDirect() {
	
	// Gewichte - Speicher
	weights.reserve(n_collo+1);
	for (int i = 0; i < n_collo+1; i++) {
		double aux = 1.;
		for (int k = 0; k < n_collo+1; k++) {
			if (i != k) {
				aux *= (time[i]-time[k]);
			}
		}
		weights.push_back( 1.0/aux );
	}
	
	// Pseudospektral Differential Matrix - Speicher
	PDM_D.reserve(n_collo);
	for(int i = 0; i < n_collo; i++) {
		PDM_D.emplace_back(n_collo+1);
	}

	#pragma omp parallel for if (n_dis > 50)
	for (int k = 0; k < n_collo; k++) {
		for (int i = 0; i < n_collo+1; i++) {
			PDM_D[k][i] = lagrangePolyDer(i,k+1);
		}
	}
}


void GaussPmTransWorhp::fillGaussPseudospectralMatrixBarycentric() {

	// Gewichte - Speicher
	weights.reserve(n_collo+1);
	for (int i = 0; i < n_collo+1; i++) {
		double aux = 1.;
		for (int k = 0; k < n_collo+1; k++) {
			if (i != k) {
				aux *= (time[i]-time[k]);
			}
		}
		weights.push_back( 1.0/aux );
	}
	
	// Pseudospektral Differential Matrix - Speicher
	PDM_D.reserve(n_collo);
	for(int i = 0; i < n_collo; i++) {
		PDM_D.emplace_back(n_collo+1);
	}
	
	for (int k = 1; k < n_collo+1; k++) {
		for (int i = 0; i < n_collo+1; i++) {
			if (i != k) {
				PDM_D[k-1][i] = (weights[i]/weights[k])/(time[k]-time[i]);
			}
		}
	}
	
	for (int k = 0; k < n_collo; k++) {
		double sum = 0.0;
		for (int j = 0; j < n_collo+1; j++) {
			if (j-1 != k) {
				sum += PDM_D[k][j];
			}
		}
		PDM_D[k][k+1] = -sum;
	}
}



double GaussPmTransWorhp::lagrangePoly(int i, double t) const {

	double aux = 1.0;
	for (int j = 1; j < n_dis-1; j++) {
		if (j != i) {
			aux *= (t - time[j])/(time[i]-time[j]);
		}
	}
	return aux;
}


double GaussPmTransWorhp::lagrangePolyDer(int i, int k) const {

	double aux = 0.0;

	for (int l = 0; l < n_dis-1; l++) {

		if (l == i) continue;

		double aux2 = 1.0;
		double aux3 = 1.0;

		for (int j = 0; j < n_dis-1; j++) {
			if (j != i && j != l) {
				aux2 *= time[k] - time[j];
			}
			if (j != i) {
				aux3 *= time[i] - time[j];
			}
			//cout << "aux " << aux2 << " " << aux3 << endl;
		}

		//cout << "aux = " << aux2/aux3 << endl;;

		aux += aux2/aux3;
	}

	return aux;
}


double GaussPmTransWorhp::interpolate(int index, double time_) const {
	
	time_ = timeToOneOne(time_);
	
	// P(t) = l(t)*sum((w_i * x(i,ode))/(t-t_i))
	// l(t) = prod(t-t_i)
	
	if (index < n_ode) { // Zustand
		
		
		double l = 1.;
		double sum = 0.;
		
		for (int i = 0; i < n_dis-1; i++) {
			
			const double y = time_ - time[i];
			if (y == 0) {
				return X[x_index(i,index)];
			}
			l *= y;
			sum += weights[i] * X[x_index(i,index)]/y;
		}
		return l*sum;
		
	} else { // Steuerung
		
		double aux = 0.0;
		for (int i = 1; i < n_dis-1; i++) {
			aux += X[u_index(i, index-n_ode)]*lagrangePoly(i, time_);
		}
		return aux;
	  
		// sollte schneller sein, klappt aber noch nicht
		/*
		double l = 1.;
		double sum = 0.;

		for (int i = 1; i < n_dis-2; i++) {
			
			const double y = time_ - time[i];
			if (y == 0) {
				return X[u_index(i,index-n_ode)];
			}
			l *= y;
			sum += weights[i]* X[u_index(i,index-n_ode)]/y;
		}

		return l*sum;
		*/
	}
}


void GaussPmTransWorhp::GetSpace(int delta1, int delta2) {

	Delta1 = delta1;
	Delta2 = delta2;

	DS_obj.useStructure(phase->obj_structure(DS_obj));
	DS_ode.useStructure(phase->ode_structure(DS_ode));
	DS_rand.useStructure(phase->rand_structure(DS_rand));
	DS_neben.useStructure(phase->neben_structure(DS_neben));

	if (n_integral) {
		
		DS_integral.useStructure(phase->integral_structure(DS_integral));
		if (DS_integral.useStructure()) {
			
			// Struktur vom Integral in Struktur der Zielfunktion uebernehmen
			for (int k = 0; k < n_integral; k++) {
				// Pruefe Abhaengigkeiten
				// Zustand
				for (int dis = 0; dis < n_dis; dis++) {
					for (int ode = 0; ode < n_ode; ode++) {
						if (DS_integral.check(k,ode)) {
							DS_obj(0, x_index(dis, ode));
						}
					}
				}
				// Steuerung
				for (int dis = 1; dis < n_dis-1; dis++) {
					for (int ctrl = 0; ctrl < n_ctrl; ctrl++) {
						if (DS_integral.check(k,n_ode+ctrl)) {
							DS_obj(0, u_index(dis, ctrl));
						}
					}
				}
				// Parameter
				for (int param = 0; param < n_param; param++) {
					if (DS_integral.check(k,n_ode + n_ctrl + param)) {
						DS_obj(0, p_index(param));
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


double GaussPmTransWorhp::x(int dis, int ode) const {
	if (dis != 0) {
		return X[n_ode + (n_ode + n_ctrl) * (dis-1) + ode];
	} else {
		return X[ode];
	}
}


double& GaussPmTransWorhp::setx(int dis, int ode) {
	if (dis != 0) {
		return X[n_ode + (n_ode + n_ctrl) * (dis-1) + ode];
	} else {
		return X[ode];
	}
}


double GaussPmTransWorhp::u(int dis, int ctrl) const {
	if (dis != 0 && dis != n_dis-1) {
		return X[n_ode + (n_ode + n_ctrl) * (dis-1) + n_ode + ctrl];
	} else {
		cout << "dont call u() with 0 or n_dis-1" << endl;
		exit(1);
		return 0.;
	}
}


double GaussPmTransWorhp::p(int param) const {
	return X[n_var - n_param + param];
}


int GaussPmTransWorhp::x_index(int dis, int ode) const {
	if (dis != 0) {
		return n_ode + (n_ode + n_ctrl) * (dis-1) + ode;
	} else {
		return ode;
	}
}


int GaussPmTransWorhp::u_index(int dis, int ctrl) const {
	if (dis != 0 && dis != n_dis-1) {
		return n_ode + (n_ode + n_ctrl) * (dis-1) + n_ode + ctrl;
	} else { // bessere Loesung?
		cout << "dont call u_index with 0 or n_dis-1" << endl;
		if (dis == 0) {
			return 1;
		} else {
			return n_dis-2;
		}
	}
}


int GaussPmTransWorhp::p_index(int param) const {
	return n_var - n_param + param;
}


void GaussPmTransWorhp::Connect(const OptVar &o, const Params &p) {

	Infty = p.Infty;

	X = &o.X[Delta1];
	X_low = &o.XL[Delta1];
	X_upp = &o.XU[Delta1];


	G = &o.G[Delta2];
	G_low = &o.GL[Delta2];
	G_upp = &o.GU[Delta2];


	Lambda = &o.Lambda[Delta1];
	Mu = &o.Mu[Delta2];

	for (int i=0; i<n_var; i++) {
		Lambda[i] = 0.0;
	}
	for (int i=0; i<n_con; i++) {
		Mu[i] = 0.0;
	}

	ZEN = &o.P[0];
	phase->zen_init(ZEN);

	std::fill(tmp_ode_1.begin(),tmp_ode_1.end(),0.0);
	std::fill(tmp_ctrl_1.begin(),tmp_ctrl_1.end(),0.0);

	for (int dis = 0; dis < n_dis; dis++) {

		phase->x_init(tmp_ode_1.data(),dis,n_dis);

		for (int ode = 0; ode < n_ode; ode++) {
			X[x_index(dis,ode)] = tmp_ode_1[ode];
		}

		if (dis != 0 && dis != n_dis-1) {
			phase->u_init(tmp_ctrl_1.data(),dis,n_dis);

			for (int ctrl = 0; ctrl < n_ctrl; ctrl++) {
				X[u_index(dis,ctrl)] = tmp_ctrl_1[ctrl];
			}
		}
	}

	phase->p_init(&X[p_index(0)]);

	phase->init();
}


void GaussPmTransWorhp::init0() {

	n_var = n_ode*(n_dis) + n_ctrl*n_collo + n_param;
	//      Kollokoation  +Endwert+  Rand  +   Nebenbed
	n_con = n_ode*n_collo + n_ode + n_rand + n_neben*n_collo;
	
	phase->localinit();

	DS_obj.Init(1, n_var);
	DS_ode.Init(n_ode, n_ode + n_ctrl + n_param);
	DS_rand.Init(n_rand, n_var);
	DS_neben.Init(n_neben, n_ode + n_ctrl + n_param);
	DS_integral.Init(n_integral, n_ode + n_ctrl + n_param);
}


void GaussPmTransWorhp::Boundary() {

	std::fill(tmp_ode_1.begin(),tmp_ode_1.end(),-Infty);
	std::fill(tmp_ode_2.begin(),tmp_ode_2.end(),+Infty);

	std::fill(tmp_ctrl_1.begin(),tmp_ctrl_1.end(),-Infty);
	std::fill(tmp_ctrl_2.begin(),tmp_ctrl_2.end(),+Infty);

	//get Boundary for state and control
	phase->x_boundary(tmp_ode_1.data(),tmp_ode_2.data());
	phase->u_boundary(tmp_ctrl_1.data(),tmp_ctrl_2.data());

	//set Boundary for state and control
	for (int dis = 0; dis < n_dis; dis++) {

		for (int ode = 0; ode < n_ode; ode++) {

			X_low[x_index(dis,ode)] = tmp_ode_1[ode];
			X_upp[x_index(dis,ode)] = tmp_ode_2[ode];
		}

		if (dis != 0 && dis != n_dis-1) {
			for (int ctrl = 0; ctrl < n_ctrl; ctrl++) {

				X_low[u_index(dis,ctrl)] = tmp_ctrl_1[ctrl];
				X_upp[u_index(dis,ctrl)] = tmp_ctrl_2[ctrl];
			}
		}
	}

	//Boundary for param
	for (int i = 0; i < n_param; i++) {
		X_low[p_index(i)] = -Infty;
		X_upp[p_index(i)] = +Infty;
	}
	phase->p_boundary( &X_low[p_index(0)], &X_upp[p_index(0)] );

	//Boundary at spec. points
	phase->var_boundary(X_low, X_upp);

	// damit die Box-Schranken durch den Startwert nicht verletzt sind
	boxConToInitGuess();

	//Boundary Collocation
	for (int i = 0; i < n_collo*n_ode+n_ode; i++) {
		G_low[i] = 0.0;
		G_upp[i] = 0.0;
	}
	
	if (n_rand) {
		std::fill(tmp_rand_1.begin(),tmp_rand_1.end(),0.0);
		std::fill(tmp_rand_2.begin(),tmp_rand_2.end(),0.0);
		
		//get Boundary for rand
		phase->rand_boundary(tmp_rand_1.data(),tmp_rand_2.data());
		
		//set Boundary for rand
		for (int i = 0; i < n_rand; i++) {
			const int index = n_collo*n_ode + n_ode + i;
			G_low[index] = tmp_rand_1[i];
			G_upp[index] = tmp_rand_2[i];
		}
	}
	
	if (n_neben) {
		std::fill(tmp_neben_1.begin(),tmp_neben_1.end(),-Infty);
		std::fill(tmp_neben_2.begin(),tmp_neben_2.end(),+Infty);
		
		//get Boundary for neben
		phase->neben_boundary(tmp_neben_1.data(),tmp_neben_2.data());
		
		//set Boundary for neben
		for (int i = 0; i < n_collo; i++) {
			for (int j = 0; j < n_neben; j++) {
				const int index = n_collo*n_ode + n_ode + n_rand + i*n_neben + j;
				G_low[index] = tmp_neben_1[j];
				G_upp[index] = tmp_neben_2[j];
			}
		}
	}
}


void GaussPmTransWorhp::Constraints2(double *GG, int /*DGflag*/) {

	// Abspeichern der ODE Auswertungen
	vector<double> tmp_ode(n_collo*n_ode);

	// Kollokationsbedingungen
	for (int dis = 1; dis < n_dis-1; dis++) {

		double *X1 = &X[x_index(dis,0)];
		double *U1 = &X[u_index(dis,0)];
		double *P  = &X[p_index(0)];

		double DX1[1000];

		phase->ode(DX1,timeToTzeroTfinal(time[dis]),X1,U1,P);
		twcount_calls.ode.f++;
		
		//ODE Wert speichern
		for (int i = 0; i < n_ode; i++) {
			tmp_ode[(dis-1)*n_ode+i] = DX1[i];
		}

		for (int ode = 0; ode < n_ode; ode++) {

			double sum = 0.0;
			//Summe
			for (int k = 0; k < n_dis-1; k++) {
				sum += PDM_D[(dis-1)][k]*x(k,ode);
			}

			// vorzeichen getauscht
			GG[(dis-1)*n_ode+ode] = 0.5*(phase->tf-phase->t0) * DX1[ode] - sum;
		}
	}
	
	//Bedingung fuer Endwert
	for (int ode = 0; ode < n_ode; ode++) {
		
		double sum = 0.0;
		for (int dis = 0; dis < n_collo; dis++) {
			sum += gauss_weights[dis]*tmp_ode[dis*n_ode + ode];
		}

		// vorzeichen getauscht
		G[n_collo*n_ode+ode] = x(n_dis-1,ode) -x(0,ode) - sum*0.5*(phase->tf-phase->t0);
	}

	// interpolation fuer Plot
	if (smoothMode) {
		// an t = -1
		for (int ode = 0; ode < n_ode; ode++) {
			displayXU[ode] = x(0,ode);
		}
		for (int ctrl = 0; ctrl < n_ctrl; ctrl++) {
			displayXU[n_ode + ctrl] = u(1,ctrl);
		}
		// an Kollokationspunkten
		for (int dis = 1; dis < displayPoints-1; dis++) {
			for (int ode = 0; ode < n_ode; ode++) {
				displayXU[dis*(n_ode+n_ctrl) + ode] = interpolate(ode,displayTime[dis]);
			}
			for (int ctrl = 0; ctrl < n_ctrl; ctrl++) {
				displayXU[dis*(n_ode+n_ctrl) + n_ode + ctrl] = interpolate(n_ode+ctrl,displayTime[dis]);
			}
		}
		//an t = +1
		for (int ode = 0; ode < n_ode; ode++) {
			displayXU[(displayPoints-1)*(n_ode+n_ctrl) + ode] = x(n_dis-1,ode);
		}
		for (int ctrl = 0; ctrl < n_ctrl; ctrl++) {
			displayXU[(displayPoints-1)*(n_ode+n_ctrl) + n_ode + ctrl] = u(n_dis-2,ctrl);
		}
	} else {
		for (int dis = 0; dis < displayPoints; dis++) {
			for (int ode = 0; ode < n_ode; ode++) {
				displayXU[dis*(n_ode+n_ctrl) + ode] = x(dis,ode);
			}
			if (dis != 0 && dis != n_dis-1) {
				for (int ctrl = 0; ctrl < n_ctrl; ctrl++) {
					displayXU[dis*(n_ode+n_ctrl) + n_ode + ctrl] = u(dis,ctrl);
				}
			} else if (dis == 0) {
				for (int ctrl = 0; ctrl < n_ctrl; ctrl++) {
					displayXU[n_ode + ctrl] = u(1,ctrl);
				}
			} else if (dis == n_dis-1) {
				for (int ctrl = 0; ctrl < n_ctrl; ctrl++) {
					displayXU[(n_dis-1)*(n_ode+n_ctrl) + n_ode + ctrl] = u(n_dis-2,ctrl);
				}
			}
		}
	}

	
	if (n_rand) {
		phase->rand( &GG[n_collo*n_ode+n_ode] );
		twcount_calls.rand.f++;
	}

	if (n_neben) {

		for (int dis = 1; dis < n_dis-1; dis++) {

			double *G1 = &GG[n_collo*n_ode+n_ode + n_rand + (dis-1)*n_neben];

			double *X1 = &X[x_index(dis,0)];
			double *U1 = &X[u_index(dis,0)];
			double *P = &X[p_index(0)];

			phase->neben(G1,timeToTzeroTfinal(time[dis]),X1,U1,P);

			twcount_calls.neben.f++;
		}
	}
}


int GaussPmTransWorhp::DF_structure(WorhpMatrix *DF, int offset) {

	DF_start = offset;
	int ind = offset;

	for (int k = 0; k < n_var; k++) {
		if (DS_obj.check(0,k)) {
			if (DF) {
				DF->row[ind] = k+1 + Delta1;
				//DF->val[ind] = 0;
			}
			ind++;
		}
	}

	DF_nnz = ind - DF_start;
	return ind;
}


void GaussPmTransWorhp::DF_calculate(WorhpMatrix &DF, double ScaleObj) {

	// Ableitung von obj zur Verfuegung?
	const bool ret = phase->obj_diff(DS_obj);

	if (ret) {
		for (int k = 0; k < DF_nnz; k++) {
			DF.val[DF_start + k] = ScaleObj * DS_obj.get(0,DF.row[DF_start+k]-1-Delta1);
		}
		twcount_calls.obj.df++;
	} else {
		cout << "keine Ableitungen von DF!" << endl;
		exit(1);
	}
}


int GaussPmTransWorhp::DG_structure(const TWfolder *f, WorhpMatrix *DG, int offset) {

	DG_start = offset;
	int ind = offset;

	for (int v = 0; v < n_var; v++) {

		//OptVar
		int var_rel;
		if (v < n_var-n_param) {
			if (v < n_ode) {
				var_rel = v;
			} else if (v > n_var-n_param-n_ode) {
				var_rel = v-(n_var-n_param-n_ode);
			} else {
				var_rel = (v-n_ode)%(n_ode+n_ctrl);
			}
		} else {
			var_rel = n_var-1-v+n_ode+n_ctrl;
		}
		
		//cout << var_rel << endl;
		
		for (int c = 0; c < n_con; c++) {

			bool add = false;

			if (c < n_ode*n_collo) { // Kollokationsbedingungen

				int con = c/n_ode; //Nebenbedingung
				int con_rel = c%n_ode; //Nebenbedingung pro Zustand

				//ODE Struktur
				if (v >= con*(n_ode+n_ctrl)+n_ode && v < (con+1)*(n_ode+n_ctrl)+n_ode) {
					if (DS_ode.check(con_rel,var_rel)) {
						add = true;
					}
				//Parameter
				} else if (v >= n_var-n_param) {
					if (DS_ode.check(con_rel,var_rel)) {
						add = true;
					}
				}

				// Summe in Kollokations-Polynom
				if (var_rel == con_rel) {
					if (v < n_var-n_param-n_ode) { //letzter Punkt kein Kollokationspunkt
						add = true;
					}
				}
			}
			else if (c < n_ode*n_collo+n_ode) { // Endbedingung
				
				int con_rel = (c-n_ode*n_collo)%n_ode; //Nebenbedingung pro Zustand
				
				if (v < n_ode) {
					if (var_rel == con_rel) {
						add = true;
					}
				} else if (v < n_var-n_param-n_ode) {
					if (DS_ode.check(con_rel,var_rel)) {
						add = true;
					}
				} else if (v < n_var-n_param) {
					//cout << v << "," << c << " : " << var_rel << " " << con_rel << endl;
					if (var_rel == con_rel) {
						add = true;
					}
				} else {
					if (DS_ode.check(con_rel,var_rel)) {
						add = true;
					}
				}
				
				//add = true;
			}
			else if (c < n_ode*n_collo+n_ode + n_rand) { // Randbedingungen

				int con_rel = (c - (n_ode*n_collo+n_ode))%n_rand;

				if (DS_rand.check(con_rel,v)) {
					add = true;
				}
				//add = true;
			}
			else if (c < n_ode*n_collo+n_ode + n_rand + n_collo*n_neben) { // Nebenbedingungen

				int con = (c - (n_ode*n_collo+n_ode + n_rand))/n_neben; //Nebenbedingung
				int con_rel = (c - (n_ode*n_collo+n_ode + n_rand))%n_neben;

				//cerr << con << " " << con_rel << endl;

				if (v >= con*(n_ode+n_ctrl)+n_ode && v < (con+1)*(n_ode+n_ctrl)+n_ode) {

					if (DS_neben.check(con_rel,var_rel)) {
						add = true;
					}

				} else if (v >= n_var-n_param) {

					if (DS_neben.check(con_rel,var_rel)) {
						add = true;
					}
				}
			}

			//add = true;
			if (add) {
				if (DG) {
					DG->col[ind] = v+1+Delta1;
					DG->row[ind] = c+1+Delta2;
					//DG->val[ind] = 0;
				}
				ind++;
			}

		}


		if(f) { //Struktur aus TWfolder
			f->DG_structure(ind, v + Delta1 + 1, DG);
		}

		//exit(1);
	}

	//exit(1);

	DG_nnz = ind - DG_start;

	return ind;
}


void GaussPmTransWorhp::DG_diff_ode(double t, int dis, int active_index) {

	DS_ode.activeindex = active_index;

	// Analytisch
	double *X1 = &X[x_index(dis, 0)];
	double *U1 = &X[u_index(dis, 0)];
	double *P  = &X[p_index(0)];

	const bool ret  = phase->ode_diff(DS_ode, t, X1, U1, P);

	if (ret) {
		twcount_calls.ode.df++;
	} else {
		cout << "keine Ableitungen von DG_ode!" << endl;
		exit(1);
	}
}


void GaussPmTransWorhp::DG_diff_ode_p(double t, int l, int dis, int active_index) {

	DS_ode.activeindex = active_index;

	// Analytisch
	double *X1 = &X[x_index(dis, 0)];
	double *U1 = &X[u_index(dis, 0)];
	double *P  = &X[p_index(0)];

	const bool ret = phase->ode_diff_p(DS_ode, t, X1, U1, P, n_ode+n_ctrl+l);

	if (ret) {
		twcount_calls.ode.dfp++;
	} else {
		cout << "keine Ableitungen von DG_ode_p!" << endl;
		exit(1);
	}
}


void GaussPmTransWorhp::DG_diff_rand() {

	const bool ret = phase->rand_diff(DS_rand);

	if (ret) {
		twcount_calls.rand.df++;
	} else {
		cout << "keine Ableitungen von DG_rand!" << endl;
		exit(1);
	}

}


void GaussPmTransWorhp::DG_diff_neben(double t, int dis) {

	if (!n_neben) return;

	// Analytisch
	double *X1 = &X[x_index(dis, 0)];
	double *U1 = &X[u_index(dis, 0)];
	double *P = &X[p_index(0)];

	const bool ret = phase->neben_diff(DS_neben, t, X1, U1, P);

	if (ret) {
		twcount_calls.neben.df++;
	} else {
		cout << "keine Ableitungen von DG_neben!" << endl;
		exit(1);
	}
}


void GaussPmTransWorhp::DG_diff_neben_p(double t, int l, int dis) {

	if (!n_neben) return;

	double *X1 = &X[x_index(dis, 0)];
	double *U1 = &X[u_index(dis, 0)];
	double *P = &X[p_index(0)];

	const bool ret = phase->neben_diff_p(DS_neben, t, X1, U1, P, n_ode+n_ctrl+l);

	if (ret) {
		twcount_calls.neben.dfp++;
	} else {
		cout << "keine Ableitungen von DG_neben_p!" << endl;
		exit(1);
	}
}


void GaussPmTransWorhp::DG_calculate(TWfolder *f, WorhpMatrix &DG) {

	// Zaehler fuer WORHP
	int index = 0;

	// damit nicht so oft ausgewertet werden muss
	int var_ode_old = -1;
	int con_ode_p_old = -1;
	int con_neben_old = -1;
	int con_neben_p_old = -1;

	if (n_rand) {
		DG_diff_rand();
	}

	for (int v = 0; v < n_var; v++) {

		//OptVar
		int var; // nur fuer PM_Matrix
		if (v < n_ode) {
			var = 0;
		} else {
			var = 1+(v-n_ode)/(n_ode+n_ctrl);
		}
		
		//cout << var << endl;
		
		int var_rel;
		if (v < n_var-n_param) {
			if (v < n_ode) {
				var_rel = v;
			} else if (v > n_var-n_param-n_ode) {
				var_rel = v-(n_var-n_param-n_ode);
			} else {
				var_rel = (v-n_ode)%(n_ode+n_ctrl);
			}
		} else {
			var_rel = n_var-1-v+n_ode+n_ctrl;
		}
		
		//cout << var << " " << var_rel << endl;
		
		for (int c = 0; c < n_con; c++) {



			if (c < n_ode*n_collo) { // Kollokationsbedingungen

				const int con = c/n_ode; //Nebenbedingung
				const int con_rel = c%n_ode; //Nebenbedingung pro Zustand

				// D[..][..] wurde schon gesetzt, fuer Diagonalen
				bool set = false;

				//ODE Struktur
				if (v >= con*(n_ode+n_ctrl)+n_ode && v < (con+1)*(n_ode+n_ctrl)+n_ode) {
					
					if (var != var_ode_old) {
						DG_diff_ode(timeToTzeroTfinal(time[var]), var, 0);
						var_ode_old = var;
					}
					
					if (DS_ode.check(con_rel, var_rel)) {
						if (var_rel == con_rel) {
							DG.val[DG_start + index] = (phase->tf-phase->t0)/2.0 * DS_ode.get(con_rel, var_rel) - PDM_D[con][var];
							set = true;
						} else {
							DG.val[DG_start + index] = (phase->tf-phase->t0)/2.0 * DS_ode.get(con_rel, var_rel);
							//cout << DS_ode.get(con_rel, var_rel) << "\t" << con_rel << "\t" << var_rel << endl;
						}
						index++;
					}
				//Parameter
				} else if (v >= n_var-n_param) {
				
					if (con != con_ode_p_old) {
						DG_diff_ode_p(timeToTzeroTfinal(time[con+1]), var_rel-(n_ode+n_ctrl), con+1, 0);
						con_ode_p_old = con;
					}
				
					if (DS_ode.check(con_rel, var_rel)) {
						DG.val[DG_start + index] = (phase->tf-phase->t0)/2.0 * DS_ode.get(con_rel, var_rel);
						index++;
					}
				}

				// Summe in Kollokations-Polynom
				if (var_rel == con_rel) {
					if (v < n_var-n_param-n_ode) { //letzter Punkt kein Kollokationspunkt
						if (!set) {
							DG.val[DG_start + index] = -1.0*PDM_D[con][var];
							index++;
						}
					}
				}
			}
			else if (c < n_ode*n_collo+n_ode) { // Endbedingung
			
				const int con_rel = (c-n_ode*n_collo)%n_ode; //Nebenbedingung pro Zustand
				
				
				if (v < n_ode) { // Startpunkt
					if (var_rel == con_rel) {
						DG.val[DG_start + index] = -1.0;
						index++;
					}
				} else if (v < n_var-n_param-n_ode) {// Zwischenpunkte
					
					if (var != var_ode_old) {
						DG_diff_ode(timeToTzeroTfinal(time[var]), var, 0);
						var_ode_old = var;
					}
					
					if (DS_ode.check(con_rel, var_rel)) {
						DG.val[DG_start + index] = -1.0*gauss_weights[var-1]*(phase->tf-phase->t0)/2.0 * DS_ode.get(con_rel, var_rel);
						index++;
					}
				} else if (v < n_var-n_param) { // Endpunkt
					//cout << v << "," << c << " : " << var_rel << " " << con_rel << endl;
					if (var_rel == con_rel) {
						DG.val[DG_start + index] = 1.0;
						index++;
					}
				} else { // freie Parameter
					
					double auxDer = 0.0;
					
					for (int i = 0; i < n_collo; i++) {
					
						DG_diff_ode_p(timeToTzeroTfinal(time[i+1]), var_rel-(n_ode+n_ctrl), i+1, 0);
						
						if (DS_ode.check(con_rel, var_rel)) {
							auxDer += gauss_weights[i]*DS_ode.get(con_rel, var_rel);
						}
					}
					
					if (DS_ode.check(con_rel, var_rel)) {
						DG.val[DG_start + index] = -1.0*(phase->tf-phase->t0)/2.0 * auxDer;
						index++;
					}
				}
				
				
			}
			else if (c < n_ode*n_collo+n_ode + n_rand) { // Randbedingungen

				const int con_rel = (c - (n_ode*n_collo+n_ode))%n_rand;
				
				if (DS_rand.check(con_rel, v)) {
					DG.val[DG_start + index] = DS_rand.get(con_rel,v);
					index++;
				}
				
			}
			else if (c < n_ode*n_collo+n_ode + n_rand + n_collo*n_neben) { // Nebenbedingungen

				const int con = (c - (n_ode*n_collo+n_ode + n_rand))/n_neben; //Nebenbedingung
				const int con_rel = (c - (n_ode*n_collo+n_ode + n_rand))%n_neben;

				//cerr << con << " " << con_rel << endl;

				if (v >= con*(n_ode+n_ctrl)+n_ode && v < (con+1)*(n_ode+n_ctrl)+n_ode) {
					
					if (con != con_neben_old) {
						DG_diff_neben(timeToTzeroTfinal(time[con+1]), con+1);
						con_neben_old = con;
					}
					
					if (DS_neben.check(con_rel, var_rel)) {
						DG.val[DG_start + index] = DS_neben.get(con_rel, var_rel);
						index++;
					}

				} else if (v >= n_var-n_param) {
				
					if (con != con_neben_p_old) {
						DG_diff_neben_p(timeToTzeroTfinal(time[con+1]), var_rel-(n_ode+n_ctrl), con+1);
						con_neben_p_old = con;
					}
					
					if (DS_neben.check(con_rel, var_rel)) {
						DG.val[DG_start + index] = DS_neben.get(con_rel, var_rel);
						index++;
					}
				}
			}
		}


		if (f) { //Struktur aus TWfolder
			f->DG_calculate(DG_start, index, v+Delta1+1, DG);
		}

		//exit(1);
	}
	//exit(1);
}


int GaussPmTransWorhp::HM_structure_ohne_Diag(int hessianstructure , WorhpMatrix& /*DF*/, WorhpMatrix& /*DG*/, WorhpMatrix* HM, int offset) {

	HM_start = offset;
	int ind = offset;

	hessianstructure = 1; //TODO
	// solange bis die (richtige) Struktur angegeben werden kann
	//cout << "hessianstructure = 1" << endl;

	if (hessianstructure == 0) { // Diagonal

	}

	else if (hessianstructure == 1) { // Full

		for (int y = 0; y < n_var; ++y) {
			for (int x = y+1; x < n_var; ++x) {
				if (HM) {
					HM->row[ind] = x + 1 + Delta1;
					HM->col[ind] = y + 1 + Delta1;
				}
				ind++;
			}
		}
	}

	else if (hessianstructure == 2) { // Blockode2, erweitert auf benachbarte Bloecke
		cout << "HM_structure 2" << endl;
		exit(1);
	}

	else if (hessianstructure == 3) { // Blockode

		/* Code aus pmTW
		const int blocksize = n_ode + n_ctrl;

		for (int y = 0; y < n_var; ++y) {

			int yblock = y / blocksize;

			for (int x = y+1; x < n_var; ++x) {

				bool use = false;

				int xblock = x / blocksize;

				if (abs(xblock-yblock) < 1) { // Id. Bloecke verbinden
					use = true;
				}

				if (x >= blocksize*n_dis ) { // mit freien Par. verbinden
					use = true;
				}

				if (use) {
					if (HM) {
						HM->row[ind] = x + 1 + Delta1;
						HM->col[ind] = y + 1 + Delta1;
					}
					ind++;
				}
			}
		}
		*/

	}

	HM_nnz_ohne_diag = ind - HM_start;

	return ind;
}


void GaussPmTransWorhp::HM_calculate1(WorhpMatrix& /*DF*/, WorhpMatrix& /*DG*/, WorhpMatrix& /*HM*/, double /*ScaleObj*/, double* /*Mu*/) {
	//hier nicht moeglich
}


void GaussPmTransWorhp::Lagrange() {

	if (!n_integral) {
		return;
	}

	std::fill_n(lagrange_integral.begin(),n_integral,0.0);

	double *P = &X[p_index(0)];
	
	for (int dis = 1; dis < n_dis-1; dis++) {

		double *X1 = &X[x_index(dis, 0)];
		double *U1 = &X[u_index(dis, 0)];

		phase->integral(tmp_integral_1.data(), timeToTzeroTfinal(time[dis]), X1, U1, P);

		for (int k = 0; k < n_integral; k++) {
			lagrange_integral[n_integral*dis + k] = lagrange_integral[n_integral*(dis-1) + k]
				+ gauss_weights[dis-1]*0.5*(phase->tf-phase->t0)*tmp_integral_1[k];
		}
	}
}


double GaussPmTransWorhp::Objective(double ScaleObj) {

	twcount_calls.obj.f++;

	Lagrange();
	double ret = phase->obj();

	for (int i = 0; i < n_integral; i++) {
		ret += lagrange_weight[i] * lagrange_integral[n_integral*n_collo+i];
	}

	return ret * ScaleObj;
}


int GaussPmTransWorhp::Integrate(int /*btableau*/) {
	return 0;
/* geht nicht(?), da Steuerung an t=-1 nicht vorhanden
	int steps = 0;

	if (btableau == 6) {

		if (butcher->stufen() == 0) {
			butcher->Init(btableau,twparameter->stepsize,false);
		}

		int startflag = 1;
		for (int i = 0; i < n_dis-1; i++) {
			int a = integrateRKF(i, startflag);
			startflag = 2;
			steps += a;

		//	cout << setw(5) << i << " " << setw(10) << a << endl;
		}
	} else {
		if (butcher->stufen() == 0) {
			butcher->Init(btableau,twparameter->stepsize);
		}

		for (int i = 0; i < n_dis-1; i++) {
			int a = integrate(i);
			steps += a;
		}
	}

#ifdef TRANSWORHP_GRAPHICS
	if (viewer) viewer->AutoScale();
#endif

	return steps;
*/
}


int GaussPmTransWorhp::integrate(int /*index*/) {
	return 0;
/*
	for (int i = 0; i < n_ode; i++) {
		tmp_ode_1[i] = x(index,i);
	}
	for (int i = 0; i < n_ctrl; i++) {
		tmp_ctrl_1[i] = u(index,i);
		tmp_ctrl_2[i] = u(index+1,i);
	}

	double h  = twparameter->stepsize;
	double t  = time[index];
	const double tend = time[index+1];

	int ct = 0;

	while (t < tend) {

		if (h+t>tend) {
			h = tend-t;
		}

		ct++;

		int ret = butcher->RungeKutta(phase,t,time[index],tend,tmp_ode_1.data(),tmp_ctrl_1.data(),tmp_ctrl_2.data(),&X[p_index(0)],h);

		if (ret == -1) {
			cout << "Too small at " << t  << endl;
			break;
		}

		if (ret == 0) {
			if (t==tend) {
				for (int l = 0; l < n_ode; l++) {
					setx(index+1,l) = tmp_ode_1[l];
				}
			}
		}
	}

	return ct;
*/
}


int GaussPmTransWorhp::integrateRKF(int /*index*/, int &/*startflag*/) {
	return 0;
/*
	double t  = time[index];
	double tend = time[index+1];

	for (int i = 0; i < n_ode; i++) {
		tmp_ode_1[i] = x(index,i);
	}
	for (int i = 0; i < n_ctrl; i++) {
		tmp_ctrl_1[i] = u(index,i);
		tmp_ctrl_2[i] = u(index+1,i);
	}

	int steps = butcher->RungeKuttaF(phase,t,time[index],tend,tmp_ode_1.data(),tmp_ctrl_1.data(),tmp_ctrl_2.data(),&X[p_index(0)],startflag);

	for (int i = 0; i < n_ode; i++){
		setx(index+1,i) = tmp_ode_1[i];
	}

	return steps;
*/
}


#ifdef TRANSWORHP_GRAPHICS
string GaussPmTransWorhp::setViewerTitle() {
	return string(" [Poly.dim = " + std::to_string(n_dis) + "] (Gauss)");
}


void GaussPmTransWorhp::updateViewer(DataStorage *ds, vector<double> &tmptime) {

	if (!smoothMode) { // diskrete Punkte
// 		ds->SetData(	n_dis, tmptime.data(),
// 				X, n_ode+n_ctrl, n_ode,
// 				&X[n_ode], n_ode+n_ctrl, n_ctrl,
// 				G,  0,  0,
// 				lagrange_integral.data(), n_integral, n_integral);
		ds->setData(	displayPoints, tmptime.data(),
				displayXU.data(), n_ode+n_ctrl, n_ode,
				displayXU.data()+n_ode, n_ode+n_ctrl, n_ctrl,
				G,  0,  0,
				lagrange_integral.data(), n_integral, n_integral);
	} else { // Polynom
		ds->setData(	displayPoints, tmptime.data(),
				displayXU.data(), n_ode+n_ctrl, n_ode,
				displayXU.data()+n_ode, n_ode+n_ctrl, n_ctrl,
				G,  0,  0,
				lagrange_integral.data(), 0, 0); // <- muss angepasst werden
	}

	if (twparameter->showGrid) {
		ds->addDisData(n_dis, n_collo, T.data(), &X[n_ode], &X[n_ode*2]);
	}
}
#endif


void GaussPmTransWorhp::PrintOCPstates(std::ostream* os) const {

	string s("OCP States");
	if (os) {
		*os << "# " << s << endl;
	} else {
		MyStatus("Viewer", s, Status::NORMAL);
	}

	s = "Phase [" + std::to_string(n_ode) + " " + std::to_string(n_ctrl) + " " + std::to_string(n_dis) + " " + std::to_string(n_param) + "] = " + std::to_string(n_var);

	if (os) {
		*os << s << endl;
	} else {
		MyStatus("OCP", s, Status::NORMAL);
	}


	if (X) {

		stringstream b;

		b << std::setw(13) << "T";
		for (int i = 0; i < n_ode; i++) {
			b << std::setw(13) << "X" << i;
		}
		for (int i = 0; i < n_ctrl; i++) {
			b << std::setw(13) << "U" << i;
		}

		if (os) {
			*os << b.str() << endl;
		} else {
			MyStatus("OCP", b.str(), Status::NORMAL);
		}

		for (int j = 0; j < n_dis; j++) {

			stringstream b;

			b <<  std::setw(13) << T[j];

			for (int i = 0; i < n_ode; i++) {
				b << std::setw(14) << x(j,i);
			}

			for (int i = 0; i < n_ctrl; i++) {
				if (j != 0 && j != n_dis-1) {
					b << std::setw(14) << u(j,i);
				} else {
					b << std::setw(14) << "-";
				}
			}

			if (os) {
				*os << b.str() << endl;
			} else {
				MyStatus("OCP", b.str(), Status::NORMAL);
			}
		}

	}


	s = "Objective function = " + std::to_string(phase->obj());

	if (os) {
		*os << s << endl;
	} else {
		MyStatus("OCP", s, Status::NORMAL);
	}
}

}
