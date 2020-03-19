#include "LobattoPmTransWORHP.h"

#include "legendreLobattoPoints.h"
#include "legendreLobattoWeights.h"

#include "TWproblem.h"
#include "TWfolder.h"
#include "butcher.h"

#include "conversion.h"

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

LobattoPmTransWorhp::LobattoPmTransWorhp(TransWorhpProblem* ph, TWparameter *twparam)
	: PmTransWorhp(ph,twparam)
{

	transworhp_type = TransWORHP_type::pseudospectral;
	
	time.resize(n_dis);
	switch (twparameter->pm_nodes_type) {
		case PMnodeType::chebyshev_lobatto:
			chebyshevLobattoTime();
			break;
		case PMnodeType::chebyshev_maxima:
			chebyshevExtremaTime();
			break;
		case PMnodeType::legendre_lobatto:
			legendreLobattoTime();
			//legendreLobattoTime2();
			break;
		default: // should never happen
			chebyshevExtremaTime();
	}
	
	/*
	for (auto ele : time) {
		cout << ele << " ";
	}
	cout << endl;
	exit(1);
	//*/
	
	for (int i = 0; i < n_dis; i++) {
		T[i] = timeToTzeroTfinal(time[i]);
	}

	// reicht hier als Groesse
	lagrange_integral.resize(n_integral*(n_dis+1));

	lobatto_weights = getLegendreLobattoWeights(n_dis);

	displayTime.reserve(displayPoints);
	for (int i = 0; i < displayPoints; i++) {
		const double t = i/(displayPoints-1.)*(phase->tf-phase->t0)+phase->t0;
		//cout << t << endl;
		displayTime.push_back(t);
	}

	displayXU.resize((n_ode+n_ctrl)*displayPoints);

	//Pseudospektral-Matrix fuellen
	fillLobattoPseudospectralMatrixBarycentric();
	//fillLobattoPseudospectralMatrixDirect();
	
	/*
	for (size_t k = 0; k < PDM_D.size(); k++) {
		for (size_t i = 0; i < PDM_D[k].size(); i++) {
			cout << std::setprecision(3) << std::setw(5) << PDM_D[k][i] << "\t";
		}
		cout << endl;
	}
	//exit(1);
	//*/
}


LobattoPmTransWorhp::~LobattoPmTransWorhp() {

}


std::string LobattoPmTransWorhp::type_G(int row) const {

	stringstream a;

	if (row < n_ode*n_dis) {
		a << " COL  " << std::setw(2) << row%n_ode << " ";
	}
	else if (row < n_ode*n_dis + n_rand) {
		a << " RAND " << std::setw(2) << (row-n_ode*n_dis)%n_rand << " ";
	}
	else if (row < n_ode*n_dis+n_rand+n_neben*n_dis) {
		a << " NB   " << std::setw(2) << (row-n_ode*n_dis-n_rand-n_neben)%n_neben << " ";
	}

	return a.str();
}


void LobattoPmTransWorhp::fromMATLAB_impl(std::ifstream &stream) {

	string line;

	/// Header
	std::getline(stream, line);
	/// Parameter p
	std::getline(stream, line);
	string line2 = string(line,2);
	vector<double> v = ToDoubleArray(line2);
	for (int i=0; i<n_param; i++) {
		X[p_index(i)] = v[i];
	}

	/// Data
	std::getline(stream, line);
	vector<double> v_last;

	v = ToDoubleArray(line);
	if (phase->freetime) v[0] = v[0]/p(0);

	v_last = v;
	for (int i = 0; i < n_dis; i++) {

		while (T[i] > v[0]) {
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

#ifdef TRANSWORHP_GRAPHICS
	if (viewer) viewer->autoScale();
#endif
}


void LobattoPmTransWorhp::chebyshevTime() {

	if (time.front() == 0) { //noch kein Wert gesetzt
		for (int i = 0; i < n_dis; i++) {
			time[i] = cos(M_PI*(2.*(n_dis-i)-1.)/(2.*n_dis));
		}
	} else {
		for (int i = 1; i < n_dis-1; i++) {
			time[i] = cos(M_PI*(2.*(n_dis-1-i)-1.)/(2.*(n_dis-2)));
		}
	}
}


void LobattoPmTransWorhp::chebyshevExtremaTime() {

	for (int i = 0; i < n_dis; i++) {
		time[i] = cos(M_PI*(n_dis-i-1.)/(n_dis-1.));
	}
}


void LobattoPmTransWorhp::chebyshevLobattoTime() {
	time.front() = -1.0;
	chebyshevTime();
	time.back() = 1.0;
}


void LobattoPmTransWorhp::legendreLobattoTime() {
	
	time = getLegendreLobattoPoints(n_dis);
	
	std::sort(time.begin(),time.end());
}


void LobattoPmTransWorhp::legendreLobattoTime2() {

	time.front() = -1.0;
	
	auto nodes = getLegendreLobattoPoints(n_dis-2);
	std::move(nodes.begin(), nodes.end(),&time[1]);
	
	time.back() = 1.0;

	std::sort(time.begin(),time.end());
}


void LobattoPmTransWorhp::fillLobattoPseudospectralMatrixDirect() {

	// Pseudospektral Differential Matrix - Speicher
	PDM_D.reserve(n_dis);
	for(int i = 0; i < n_dis; i++) {
		PDM_D.emplace_back(n_dis);
	}

	// Gewichte - Speicher
	weights.reserve(n_dis);
	for (int i = 0; i < n_dis; i++) {
		double aux = 1.;
		for (int k = 0; k < n_dis; k++) {
			if (i != k) {
				aux *= (time[i]-time[k]);
			}
		}
		weights.push_back( 1.0/( aux ) );
	}

	#pragma omp parallel for if (n_dis > 50)
	for (int k = 0; k < n_dis; k++) {
		for (int i = 0; i < n_dis; i++) {
			PDM_D[k][i] = lagrangePolyDer(i,k);
		}
	}
}


void LobattoPmTransWorhp::fillLobattoPseudospectralMatrixBarycentric() {

	// Gewichte - Speicher
	weights.reserve(n_dis);
	for (int i = 0; i < n_dis; i++) {
		double aux = 1.;
		for (int j = 0; j < n_dis; j++) {
			if (i != j) {
				aux *= (time[i]-time[j]);
			}
		}
		weights.push_back( 1.0/( aux ) );
	}

	// Pseudospektral Differential Matrix - Speicher
	PDM_D.reserve(n_dis);
	for(int i = 0; i < n_dis; i++) {
		PDM_D.emplace_back(n_dis);
	}

	for (int i = 0; i < n_dis; i++) {
		for (int k = 0; k < n_dis; k++) {
			if (i != k) {
				PDM_D[k][i] = (weights[i]/weights[k])/(time[k]-time[i]);
			}
		}
	}
	
	for (int k = 0; k < n_dis; k++) {
		double sum = 0.0;
		for (int j = 0; j < n_dis; j++) {
			if (k != j) {
				sum += PDM_D[k][j];
			}
		}
		PDM_D[k][k] = -sum;
	}
}


double LobattoPmTransWorhp::lagrangePoly(int i, double t) const {

	double aux = 1.0;
	for (int j = 0; j < n_dis; j++) {
		if (j != i) {
			aux *= (t - time[j])/(time[i]-time[j]);
		}
	}
	return aux;
}


double LobattoPmTransWorhp::lagrangePolyDer(int i, int k) const {

	double aux = 0.0;

	for (int l = 0; l < n_dis; l++) {

		if (l == i) continue;

		double aux2 = 1.0;
		double aux3 = 1.0;

		for (int j = 0; j < n_dis; j++) {
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


double LobattoPmTransWorhp::interpolate(int index, double time_) const {
	
	time_ = timeToOneOne(time_);
	
	// P(t) = l(t)*sum((w_i * x(i,ode))/(t-t_i))
	// l(t) = prod(t-t_i)
	
	double l = 1.;
	double sum = 0.;
	
	for (int i = 0; i < n_dis; i++) {
		
		const double y = time_ - time[i];
		if (y == 0) {
			return X[(n_ode + n_ctrl) * i + index];
		}
		l *= y;
		sum += weights[i]* X[(n_ode + n_ctrl) * i + index]/y;
	}
	return l*sum;
	
	
	// direkte Formel
	/*
	double aux = 0.0;
	for (int i = 0; i < n_dis; i++) {
		aux += X[(n_ode + n_ctrl) * i + index]*lagrangePoly(i,time_);
	}
	return aux;
	*/
}


void LobattoPmTransWorhp::GetSpace(int delta1, int delta2) {

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
				// Zustand + Steuerung
				for (int dis = 0; dis < n_dis; dis++) {
					for (int ode = 0; ode < n_ode; ode++) {
						if (DS_integral.check(k,ode)) {
							DS_obj(0, x_index(dis, ode));
						}
					}
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

double LobattoPmTransWorhp::x(int dis, int ode) const {
	return X[(n_ode + n_ctrl) * dis + ode];
}


double &LobattoPmTransWorhp::setx(int dis, int ode) {
	return X[(n_ode + n_ctrl) * dis + ode];
}


double LobattoPmTransWorhp::u(int dis, int ctrl) const {
	return X[(n_ode + n_ctrl) * dis + n_ode + ctrl];
}


double LobattoPmTransWorhp::p(int param) const {
	return X[(n_ode + n_ctrl) * n_dis + param];
}


int LobattoPmTransWorhp::x_index(int dis, int ode) const {
	return (n_ode + n_ctrl) * dis + ode;
}


int LobattoPmTransWorhp::u_index(int dis, int ctrl) const {
	return (n_ode + n_ctrl) * dis + n_ode + ctrl;
}


int LobattoPmTransWorhp::p_index(int param) const {
	return (n_ode + n_ctrl) * n_dis + param;
}


void LobattoPmTransWorhp::Connect(const OptVar &o, const Params &p) {

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

		const int index1 = (n_ode + n_ctrl) * dis;
		for (int ode = 0; ode < n_ode; ode++) {
			X[index1+ode] = tmp_ode_1[ode];
		}

		phase->u_init(tmp_ctrl_1.data(),dis,n_dis);

		const int index2 = (n_ode + n_ctrl) * dis + n_ode;
		for (int ctrl = 0; ctrl < n_ctrl; ctrl++) {
			X[index2+ctrl] = tmp_ctrl_1[ctrl];
		}
	}

	phase->p_init(&X[n_dis*(n_ctrl+n_ode)]);

	phase->init();
}


void LobattoPmTransWorhp::init0() {

	n_var = (n_ode + n_ctrl) * n_dis + n_param;
	n_con = n_ode * n_dis + n_rand + n_neben*n_dis;

	phase->localinit();

	DS_obj.Init(1, n_var);
	DS_ode.Init(n_ode, n_ode + n_ctrl + n_param);
	DS_rand.Init(n_rand, n_var);
	DS_neben.Init(n_neben, n_ode + n_ctrl + n_param);
	DS_integral.Init(n_integral, n_ode + n_ctrl + n_param);
}


void LobattoPmTransWorhp::Boundary() {

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

			const int index = (n_ode + n_ctrl) * dis + ode;
			X_low[index] = tmp_ode_1[ode];
			X_upp[index] = tmp_ode_2[ode];
		}

		for (int ctrl = 0; ctrl < n_ctrl; ctrl++) {

			const int index = (n_ode + n_ctrl) * dis + n_ode + ctrl;
			X_low[index] = tmp_ctrl_1[ctrl];
			X_upp[index] = tmp_ctrl_2[ctrl];
		}
	}

	//Boundary for param
	for (int i = 0; i < n_param; i++) {
		X_low[n_dis*(n_ctrl+n_ode)+i] = -Infty;
		X_upp[n_dis*(n_ctrl+n_ode)+i] = +Infty;
	}
	phase->p_boundary( &X_low[n_dis*(n_ctrl+n_ode)], &X_upp[n_dis*(n_ctrl+n_ode)] );

	//Boundary at spec. points
	phase->var_boundary(X_low, X_upp);

	// damit die Box-Schranken durch den Startwert nicht verletzt sind
	boxConToInitGuess();

	for (int i = 0; i < n_dis*n_ode; i++) {
		G_low[i] = 0.0;
		G_upp[i] = 0.0;
	}

	std::fill(tmp_rand_1.begin(),tmp_rand_1.end(),0.0);
	std::fill(tmp_rand_2.begin(),tmp_rand_2.end(),0.0);

	//get Boundary for rand
	if (n_rand) phase->rand_boundary(tmp_rand_1.data(),tmp_rand_2.data());

	//set Boundary for rand
	for (int i = 0; i < n_rand; i++) {
		const int index = n_dis * n_ode + i;
		G_low[index] = tmp_rand_1[i];
		G_upp[index] = tmp_rand_2[i];
	}

	std::fill(tmp_neben_1.begin(),tmp_neben_1.end(),-Infty);
	std::fill(tmp_neben_2.begin(),tmp_neben_2.end(),+Infty);

	//get Boundary for neben
	if (n_neben) phase->neben_boundary(tmp_neben_1.data(),tmp_neben_2.data());

	//set Boundary for neben
	for (int i = 0; i < n_dis; i++) {
		for (int j = 0; j < n_neben; j++) {
			const int index = n_dis * n_ode + n_rand + i*n_neben + j;
			G_low[index] = tmp_neben_1[j];
			G_upp[index] = tmp_neben_2[j];
		}
	}
}


void LobattoPmTransWorhp::Constraints2(double *GG, int /*DGflag*/) {

	for (int dis = 0; dis < n_dis; dis++) {

		double *X1 = &X[(n_ode+n_ctrl)*dis];
		double *U1 = &X[(n_ode+n_ctrl)*dis + n_ode];
		double *P  = &X[(n_ode+n_ctrl)*n_dis];

		double DX1[1000];

		phase->ode(DX1,timeToTzeroTfinal(time[dis]),X1,U1,P);
		twcount_calls.ode.f++;

		for (int ode = 0; ode < n_ode; ode++) {

			double sum = 0.0;

			//Summe
			for (int k = 0; k < n_dis; k++) {
				sum += PDM_D[dis][k]*x(k,ode);
			}

			GG[dis*n_ode+ode] = 0.5*(phase->tf-phase->t0) * DX1[ode] - sum;
		}
	}

	if (smoothMode) {
		for (int dis = 0; dis < displayPoints; dis++) {

			for (int ode = 0; ode < n_ode; ode++) {
				displayXU[dis*(n_ode+n_ctrl) + ode] = interpolate(ode,displayTime[dis]);
			}
			for (int ctrl = 0; ctrl < n_ctrl; ctrl++) {
				displayXU[dis*(n_ode+n_ctrl) + n_ode + ctrl] = interpolate(n_ode+ctrl,displayTime[dis]);
			}
		}
	}

	if (n_rand) {
		phase->rand( &GG[n_dis*n_ode] );
		twcount_calls.rand.f++;
	}

	if (n_neben) {

		for (int dis = 0; dis < n_dis; dis++) {

			double *G1 = &GG[n_dis*n_ode + n_rand + dis*n_neben];

			double *X1 = &X[(n_ode+n_ctrl)*dis];
			double *U1 = &X[(n_ode+n_ctrl)*dis + n_ode];
			double *P = &X[(n_ode+n_ctrl)*n_dis];

			phase->neben(G1,timeToTzeroTfinal(time[dis]),X1,U1,P);

			twcount_calls.neben.f++;
		}
	}
}


int LobattoPmTransWorhp::DF_structure(WorhpMatrix *DF, int offset) {

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


void LobattoPmTransWorhp::DF_calculate(WorhpMatrix &DF, double ScaleObj) {

	// Ableitung von obj zur Verfuegung?
	bool ret = phase->obj_diff(DS_obj);

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


int LobattoPmTransWorhp::DG_structure(const TWfolder *f, WorhpMatrix *DG, int offset) {

	DG_start = offset;
	int ind = offset;

	for (int v = 0; v < n_var; v++) {

		//cout << var << endl;

		int var_rel; //OptVar
		if (v < n_var-n_param) {
			var_rel = v%(n_ode+n_ctrl);
		} else {
			var_rel = v-(n_ode+n_ctrl)*n_dis + n_ode+n_ctrl;
		}

		for (int c = 0; c < n_con; c++) {

			bool add = false;

			if (c < n_ode*n_dis) { // Kollokationsbedingungen

				int con = c/n_ode; //Nebenbedingung
				int con_rel = c%n_ode; //Nebenbedingung pro Zustand

				//cerr << con_rel << " " << var_rel << " " << con << " " << "" << endl;
				//cerr << v << " " << con*(n_ode+n_ctrl) << endl;
				if (v >= con*(n_ode+n_ctrl) && v < (con+1)*(n_ode+n_ctrl)) {

					if (DS_ode.check(con_rel,var_rel)) {
						add = true;
					}

				} else if (v >= n_var-n_param) {

					if (DS_ode.check(con_rel,var_rel)) {
						add = true;
					}
				}

				// Summe in Kollokations-Polynom
				if (var_rel == con_rel) {
					add = true;
				}
			}
			else if (c < n_ode*n_dis + n_rand) { // Randbedingungen

				int con_rel = (c - n_dis*n_ode)%n_rand;

				if (DS_rand.check(con_rel,v)) {
					add = true;
				}

			}
			else if (c < n_ode*n_dis + n_rand + n_dis*n_neben) { // Nebenbedingungen

				int con = (c - n_dis*n_ode - n_rand)/n_neben; //Nebenbedingung
				int con_rel = (c - n_dis*n_ode - n_rand)%n_neben;

				//cerr << con << " " << con_rel << endl;

				if (v >= con*(n_ode+n_ctrl) && v < (con+1)*(n_ode+n_ctrl)) {

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


void LobattoPmTransWorhp::DG_diff_ode(double t, int dis, int active_index) {

	DS_ode.activeindex = active_index;

	// Analytisch
	double *X1 = &X[(n_ode+n_ctrl)* dis];
	double *U1 = &X[(n_ode+n_ctrl)* dis + n_ode];
	double *P = &X[(n_ode+n_ctrl)*n_dis];

	bool ret  = phase->ode_diff(DS_ode, t, X1, U1, P);
	if (ret) {
		twcount_calls.ode.df++;
	} else {
		cout << "keine Ableitungen von DG_ode!" << endl;
		exit(1);
	}

}


void LobattoPmTransWorhp::DG_diff_ode_p(double t, int l, int dis, int active_index) {

	DS_ode.activeindex = active_index;

	// Analytisch
	double *X1 = &X[(n_ode+n_ctrl)*dis];
	double *U1 = &X[(n_ode+n_ctrl)*dis + n_ode];
	double *P = &X[(n_ode+n_ctrl)*n_dis];

	bool ret = phase->ode_diff_p(DS_ode, t, X1, U1, P, n_ode+n_ctrl+l);
	if (ret) {
		twcount_calls.ode.dfp++;
	} else {
		cout << "keine Ableitungen von DG_ode_p!" << endl;
		exit(1);
	}
}


void LobattoPmTransWorhp::DG_diff_rand() {

	bool ret = phase->rand_diff(DS_rand);

	if (ret) {
		twcount_calls.rand.df++;
	} else {
		cout << "keine Ableitungen von DG_rand!" << endl;
		exit(1);
	}

}


void LobattoPmTransWorhp::DG_diff_neben(double t, int dis) {

	if (!n_neben) return;

	// Analytisch
	double *X1 = &X[(n_ode+n_ctrl)*dis];
	double *U1 = &X[(n_ode+n_ctrl)*dis + n_ode];
	double *P = &X[(n_ode+n_ctrl)*n_dis];

	bool ret = phase->neben_diff(DS_neben, t, X1, U1, P);
	if (ret) {
		twcount_calls.neben.df++;
	} else {
		cout << "keine Ableitungen von DG_neben!" << endl;
		exit(1);
	}
}


void LobattoPmTransWorhp::DG_diff_neben_p(double t, int l, int dis) {

	if (!n_neben) return;

	double *X1 = &X[(n_ode+n_ctrl)*dis];
	double *U1 = &X[(n_ode+n_ctrl)*dis + n_ode];
	double *P = &X[(n_ode+n_ctrl)*n_dis];

	bool ret = phase->neben_diff_p(DS_neben, t, X1, U1, P, n_ode+n_ctrl+l);
	if (ret) {
		twcount_calls.neben.dfp++;
	} else {
		cout << "keine Ableitungen von DG_neben_p!" << endl;
		exit(1);
	}
}


void LobattoPmTransWorhp::DG_calculate(TWfolder *f, WorhpMatrix &DG) {

	// Zaehler fuer WORHP
	int index = 0;

	// damit nicht so oft ausgewertet werden muss
	int con_ode_old = -1;
	int con_ode_p_old = -1;
	int con_neben_old = -1;
	int con_neben_p_old = -1;

	if (n_rand) {
		DG_diff_rand();
	}

	for (int v = 0; v < n_var; v++) {

		const int var = v/(n_ode+n_ctrl); // nur fuer PM_Matrix
		int var_rel; //OptVar
		if (v < n_var-n_param) {
			var_rel = v%(n_ode+n_ctrl);
		} else {
			var_rel = v-(n_ode+n_ctrl)*n_dis + n_ode+n_ctrl;
		}

		for (int c = 0; c < n_con; c++) {


			if (c < n_ode*n_dis) { // Kollokationsbedingungen


				int con = c/n_ode; //Nebenbedingung
				int con_rel = c%n_ode; //Nebenbedingung pro Zustand

				//cerr << con_rel << " " << var_rel << " " << con << " " << "" << endl;
				//cerr << v << " " << con*(n_ode+n_ctrl) << endl;

				// D[..][..] wurde schon gesetzt, fuer Diagonalen
				bool set = false;

				if (v >= con*(n_ode+n_ctrl) && v < (con+1)*(n_ode+n_ctrl)) {

					if (con != con_ode_old) {
						DG_diff_ode(timeToTzeroTfinal(time[con]),con,0);
						con_ode_old = con;
					}

					//cerr << con << endl;

					if (DS_ode.check(con_rel,var_rel)) {

						if (var_rel == con_rel) {
							DG.val[DG_start + index] = (phase->tf-phase->t0)/2.0 * DS_ode.get(con_rel,var_rel) - PDM_D[con][var];
							set = true;
						} else {
							DG.val[DG_start + index] = (phase->tf-phase->t0)/2.0 * DS_ode.get(con_rel,var_rel);
						}
						index++;
					}

				} else if (v >= n_var-n_param) {

					if (con != con_ode_p_old) {
						DG_diff_ode_p(timeToTzeroTfinal(time[con]),var_rel-(n_ode+n_ctrl),con,0);
						con_ode_p_old = con;
					}

					if (DS_ode.check(con_rel,var_rel)) {
						DG.val[DG_start + index] = (phase->tf-phase->t0)/2.0 * DS_ode.get(con_rel,var_rel);
						index++;
					}
				}

				// Summe in Kollokations-Polynom
				if (var_rel == con_rel) {
					if (!set) {
						DG.val[DG_start + index] = -1.0*PDM_D[con][var];
						index++;
					}
				}
			}
			else if (c < n_ode*n_dis + n_rand) { // Randbedingungen

				int con_rel = (c - n_dis*n_ode)%n_rand;

				if (DS_rand.check(con_rel,v)) {
					DG.val[DG_start + index] = DS_rand.get(con_rel,v);
					index++;
				}

			}
			else if (c < n_ode*n_dis + n_rand + n_dis*n_neben) { // Nebenbedingungen

				int con = (c - n_dis*n_ode - n_rand)/n_neben; //Nebenbedingung
				int con_rel = (c - n_dis*n_ode - n_rand)%n_neben;

				//cerr << con << " " << con_rel << endl;
				//cerr << con << endl;

				if (v >= con*(n_ode+n_ctrl) && v < (con+1)*(n_ode+n_ctrl)) {

					if (con != con_neben_old) {
						DG_diff_neben(timeToTzeroTfinal(time[con]),con);
						con_neben_old = con;
					}

					if (DS_neben.check(con_rel,var_rel)) {
						DG.val[DG_start + index] = DS_neben.get(con_rel,var_rel);
						index++;
					}

				} else if (v >= n_var-n_param) {

					if (con != con_neben_p_old) {
						DG_diff_neben_p(timeToTzeroTfinal(time[con]),var_rel-(n_ode+n_ctrl),con);
						con_neben_p_old = con;
					}

					if (DS_neben.check(con_rel,var_rel)) {
						DG.val[DG_start + index] = DS_neben.get(con_rel,var_rel);
						index++;
					}
				}
			}
		}

		if (f) {
			f->DG_calculate(DG_start,index,v+Delta1+1,DG);
		}

		//exit(1);
	}
	//exit(1);

}


int LobattoPmTransWorhp::HM_structure_ohne_Diag(int hessianstructure , WorhpMatrix& /*DF*/, WorhpMatrix& /*DG*/, WorhpMatrix* HM, int offset) {

	HM_start = offset;
	int ind = offset;

	//hessianstructure = 1; //TODO
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

	}

	HM_nnz_ohne_diag = ind - HM_start;

	return ind;
}


void LobattoPmTransWorhp::HM_calculate1(WorhpMatrix& /*DF*/, WorhpMatrix& /*DG*/, WorhpMatrix& /*HM*/, double /*ScaleObj*/, double* /*Mu*/) {
	//hier nicht moeglich
}


void LobattoPmTransWorhp::Lagrange() {

	if (!n_integral) {
		return;
	}

	std::fill_n(lagrange_integral.begin(),n_integral,0.0);

	double *P = &X[p_index(0)];
	
	for (int dis = 0; dis < n_dis; dis++) {

		double *X1 = &X[x_index(dis, 0)];
		double *U1 = &X[u_index(dis, 0)];

		phase->integral(tmp_integral_1.data(), timeToTzeroTfinal(time[dis]), X1, U1, P);

		for (int k = 0; k < n_integral; k++) {
			lagrange_integral[n_integral*(dis+1) + k] = lagrange_integral[n_integral*dis + k] + lobatto_weights[dis]*0.5*(phase->tf-phase->t0)*tmp_integral_1[k];
		}
	}
}


double LobattoPmTransWorhp::Objective(double ScaleObj) {

	twcount_calls.obj.f++;

	Lagrange();
	double ret = phase->obj();

	for (int i = 0; i < n_integral; i++) {
		ret += lagrange_weight[i] * lagrange_integral[n_integral*n_dis+i];
	}

	return ret * ScaleObj;
}


int LobattoPmTransWorhp::Integrate(int btableau) {

	int steps = 0;

	if (btableau==6) {

		if (butcher.stufen() == 0) {
			butcher.Init(btableau, twparameter->stepsize, false);
		}

		int startflag = 1;
		for (int i = 0; i < n_dis-1; i++) {
			const int a = integrateRKF(i, startflag);
			startflag = 2;
			steps += a;

		//	cout << setw(5) << i << " " << setw(10) << a << endl;
		}

	}
	else {

		if (butcher.stufen() == 0) {
			butcher.Init(btableau, twparameter->stepsize);
		}

		for (int i = 0; i < n_dis-1; i++) {
			int a = integrate(i);
			steps += a;
		}
	}

#ifdef TRANSWORHP_GRAPHICS
	if (viewer) viewer->autoScale();
#endif

	return steps;

}


int LobattoPmTransWorhp::integrate(int index) {

	for (int l=0; l<n_ode; l++) {
		tmp_ode_1[l] = x(index,l);
	}
	for (int l=0; l<n_ctrl; l++) {
		tmp_ctrl_1[l] = u(index,l);
		tmp_ctrl_2[l] = u(index+1,l);
	}

	double h  = twparameter->stepsize;
	double t  = timeToTzeroTfinal(time[index]);
	const double tend = timeToTzeroTfinal(time[index+1]);

	int ct = 0;

	while (t < tend) {

		if (h+t>tend) {
			h = tend-t;
		}

		ct++;

		int ret = butcher.RungeKutta(phase,t,timeToTzeroTfinal(time[index]),tend,tmp_ode_1.data(),tmp_ctrl_1.data(),tmp_ctrl_2.data(),&X[p_index(0)],h);

		if (ret==-1) {
			cout << "Too small at " << t  << endl;
			break;
		}

		if (ret==0) {
			if (t==tend) {
				for (int l = 0; l < n_ode; l++) {
					setx(index+1,l) = tmp_ode_1[l];
				}
			}
		}
	}

	return ct;
}


int LobattoPmTransWorhp::integrateRKF(int index, int &startflag) {

	double t  = timeToTzeroTfinal(time[index]);
	double tend = timeToTzeroTfinal(time[index+1]);

	for (int l=0; l<n_ode; l++) {
		tmp_ode_1[l] = x(index,l);
	}
	for (int l=0; l<n_ctrl; l++) {
		tmp_ctrl_1[l] = u(index,l);
		tmp_ctrl_2[l] = u(index+1,l);
	}

	int steps = butcher.RungeKuttaF(phase,t,timeToTzeroTfinal(time[index]),tend,tmp_ode_1.data(),tmp_ctrl_1.data(),tmp_ctrl_2.data(),&X[p_index(0)],startflag);

	for (int l=0; l<n_ode; l++){
		setx(index+1,l) = tmp_ode_1[l];
	}

	return steps;
}


#ifdef TRANSWORHP_GRAPHICS
string LobattoPmTransWorhp::setViewerTitle() {
	return string(" [Poly.dim = " + std::to_string(n_dis) + "] (Lobatto)");
}


void LobattoPmTransWorhp::updateViewer(DataStorage *ds, vector<double> &tmptime) {

	if (!smoothMode) { // diskrete Punkte
		ds->setData(	n_dis, tmptime.data(),
				X, n_ode+n_ctrl, n_ode,
				&X[n_ode], n_ode+n_ctrl, n_ctrl,
				G,  0,  0,
				lagrange_integral.data(), n_integral, n_integral);
	} else { // Polynom
		ds->setData(	displayPoints, tmptime.data(),
				displayXU.data(), n_ode+n_ctrl, n_ode,
				displayXU.data()+n_ode, n_ode+n_ctrl, n_ctrl,
				G,  0,  0,
				lagrange_integral.data(), n_integral, n_integral);
	}

	if (twparameter->showGrid) {
		ds->addDisData(n_dis, n_dis, T.data(), X, &X[n_ode]);
	}
}
#endif

}
