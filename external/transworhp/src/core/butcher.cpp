#include "butcher.h"

#include "butcherTableau.h"

#include "worhp/worhp.h"

#include "TransWORHP.h"
#include "TWparameter.h"
#include "TWproblem.h"
#include "twstatus.h"

#include "conversion.h"

#include "rkf45.hpp"

#include <iostream>
#include <iomanip>
#include <fstream>

#ifdef TW_WITH_BOOST
#include <boost/numeric/odeint.hpp>
using namespace boost::numeric::odeint;
#endif

namespace tw {

#if defined(TW_WITH_SUPERLU) || defined(WITH_LAPACK)
#include "ExplTransWORHP.h"
#endif

#ifdef TW_WITH_SUPERLU
extern "C"
{
#include "superlu/slu_ddefs.h"
}
#endif


#ifdef WITH_LAPACK
extern "C"
{
	void dgesv_(int &n, int &nrhs, double *A, int &lda, int *ipiv, double *b, int &ldb, int &info);
	void dgetrs_(char *trans, int &n, int &nrhs, double *A,int &lda, int *ipiv, double *b, int &ldb, int &info);
	void dgetrf_(int &m, int &n, double* A, int &lda, int *ipiv, int &info);
}
#endif


using namespace std;


Butcher::Butcher(int n_ode, int n_ctrl, double abserr, double relerr) :
	h0(.0),
	isInit(false),
	amin(.5), amax(1.5), aa(0.9),
	maxEpsEst(1.0), minEpsEst(1.0), eps(1e-2),
	hmin(1e-10),
	abserr(abserr), relerr(relerr),
	tmp_x(n_ode), tmp_u(n_ctrl),
	tmp_ret1(n_ode), tmp_ret2(n_ode),
	s(0), p(0),
	n_ode(n_ode), n_ctrl(n_ctrl) {


#if defined(TW_WITH_SUPERLU) || defined(WITH_LAPACK)
	// Struktur fuer ROW-Verfahren
	ROW_DS.useStructure(true);
	ROW_DS.finish();
	ROW_DS.Init(n_ode,n_ode);

	ROW_DS_check = false;
#endif
}

void Butcher::Init(int index, double stepsize, bool verbose) {

	isInit = true;

#if defined(TW_WITH_SUPERLU) || defined(WITH_LAPACK)
	// ROW-Verfahren
	if (index > 99) {
		s = 2;
		p = 3;

#ifdef TW_WITH_SUPERLU
		MyStatus("Butcher", "ROW-Method (SuperLU)", Status::WARN);
#else
		MyStatus("Butcher", "ROW-Method (Lapack)", Status::WARN);
#endif

		// Butcher-Tableau fuer ROW-Verfahren (Wiki)
		// https://de.wikipedia.org/wiki/Rosenbrock-Wanner-Verfahren
		///////////////////////////////////
		a[0][0] = 0;
		a[0][1] = 0;
		a[1][0] = 2.0/3.0;
		a[1][1] = 0;

		b[0] = 1.0/4.0;
		b[1] = 3.0/4.0;

		c[0] = 0;
		c[1] = 2.0/3.0;

		const double g = (1.0+1.0/sqrt(3.0))/2.0;

		gamma[0][0] = g;
		gamma[0][1] = 0;
		gamma[1][0] = -4.0/3.0*g;
		gamma[1][1] = g;

		d[0] = gamma[0][0];
		d[1] = gamma[0][1]+gamma[1][1];
		///////////////////////////////////

		nnzOfA = 0;

		return;
	}
#endif


	//boost integrate
	if (index >= 20 && index < 30) {
		MyStatus("Butcher", "Tableau " + names[index], Status::WARN);

#ifdef TW_WITH_BOOST
		h0 = stepsize;
		return;
#else
		MyStatus("Butcher", "Please activate boost.", Status::ERR);
#endif
	}

	if (static_cast<size_t>(index) >= tableaus.size()) {
		MyStatus("Butcher", "Butcher-Tableau ID is not valid!", Status::ERR);
		MyStatus("Butcher", "chaging to Euler-Method!", Status::ERR);
		index = 10;
	}

	h0 = stepsize;
	stringstream is(tableaus[index]);

	if (is) {
		is >> *this;

		string str;
		if (hasStepSize()==0) str = "  (" + to_string(h0) + ")";
		MyStatus("Butcher", "Tableau " + names[index] + std::move(str), Status::WARN);

		if (verbose) {
			stringstream bb;
			bb << *this;
			MyStatus("Butcher", bb.str(), Status::NORMAL);
		}
	} else {
		MyStatus("Butcher", "Tableau " + to_string(index) + " not valid", Status::ERR);
	}

	if (hasStepSize()) {
		for (int i = 0; i < p; i++) {
			maxEpsEst *= amax/aa;
			minEpsEst *= amin/aa;
		}
	}

	k.reserve(s);
	for (int i = 0; i < s; i++) {
		k.emplace_back(n_ode);
	}
}

double Butcher::getratio(istream &os) {

	string s;
	os >> s;
	vector<int> b = ToIntArray(s,"/");
	if (b.size() == 2) {
		return (b[0]/static_cast<double>(b[1]));
	} else {
		return ToDouble(s);
	}
}


istream &operator>>(istream &os, class Butcher &b) {

	os >> b.s;
	os >> b.p;
	for (int i=0; i<b.s; i++) {
		b.c[i] = b.getratio(os);

		for (int j=0; j<i; j++) {
			b.a[i][j] = b.getratio(os);
		}
	}
	for (int i=0; i<b.s; i++) {
		b.b[i] = b.getratio(os);
	}

	if (b.hasStepSize()) {
		for (int i=0; i<b.s; i++) {
			b.b2[i] = b.getratio(os);
		}
	}

	return os;
}

ostream &operator<<(ostream &os, const Butcher &b) {

	os.setf(ios::fixed);
	os << setprecision(7);

	for (int i=0; i<b.s; i++) {
		os << setw(14) << b.c[i] << " | ";
		for (int j=0; j<i; j++) {
			os << setw(14)<< b.a[i][j];
		}
		os << endl;
	}
	os << "---------------+-";
	for (int i=0; i<b.s; i++) {
		os << "--------------";
	}
	os << endl;
	os << "               | " ;
	for (int i=0; i<b.s; i++) {
		os << setw(14) << b.b[i];
	}
	os << endl;
	if (b.hasStepSize()) {
		os << "               | " ;
		for (int i=0; i<b.s; i++) {
			os << setw(14) << b.b2[i];
		}
		os << endl;
	}
	os << endl;

	return os;
}

TransWorhpProblem *ph=nullptr;

void f ( double t, double y[], double yp[], double t0, double tend, double *u1, double *u2, double *param, int n_ctrl) {

	double tmp_u[100];

	double pos = (t-t0)/(tend-t0);

	for (int l=0; l<n_ctrl; l++) {
		if (ph->solver->twparameter->linInter==1) { // lineare Interpolation der Steuerung
			tmp_u[l] = u1[l] + (u2[l]-u1[l]) * pos;
		} else { // konstante Steuerung
			tmp_u[l] = u1[l];
		}
	}

	ph->ode(yp,t,y, tmp_u,param);
}


#ifdef TW_WITH_BOOST
template<class stepper_type>
int Butcher::integrate_with_boost(TransWorhpProblem* tw, double t0, double tend, std::vector<double> &x, const double* param) const {

	stepper_type stepper;

	const size_t steps = integrate_const(stepper,
		[&tw, &param](const std::vector<double> &x , std::vector<double> &dxdt , double t) {

			std::vector<double> tmp_u(tw->n_ctrl);
			tw->solver->GetControl(tmp_u.data(), t);

			tw->ode(dxdt.data(), t, x.data(), tmp_u.data(), param);

		}, x, t0, tend, h0/*,
		[](const std::vector<double> &x , double t) {
			std::cout << t << std::endl;
		}*/);

	return static_cast<int>(steps);
}


template<class stepper_type>
int Butcher::integrate_with_boost_controlled(TransWorhpProblem* tw, double t0, double tend, std::vector<double> &x, const double* param) const {

	stepper_type stepper;
	auto controlled_stepper = make_controlled(abserr, relerr, stepper);

	const size_t steps = integrate_const(controlled_stepper,
		[&tw, &param](const std::vector<double> &x , std::vector<double> &dxdt , double t) {

			std::vector<double> tmp_u(tw->n_ctrl);
			tw->solver->GetControl(tmp_u.data(), t);

			tw->ode(dxdt.data(), t, x.data(), tmp_u.data(), param);

		}, x, t0, tend, h0/*,
		[](const std::vector<double> &x , double t) {
			std::cout << t << std::endl;
		}*/);

	return static_cast<int>(steps);
}


int Butcher::integrate_with_boost(TransWorhpProblem* tw, double t0, double tend, std::vector<double> &x, const double* param, int butchertableau) const {

	const size_t multistep_order = 6;

	switch (butchertableau) {
		case 20:
			return integrate_with_boost<euler<std::vector<double>>>(tw, t0, tend, x, param);
			break;
		case 21:
			return integrate_with_boost<modified_midpoint<std::vector<double>>>(tw, t0, tend, x, param);
			break;
		case 22:
			return integrate_with_boost<runge_kutta4<std::vector<double>>>(tw, t0, tend, x, param);
			break;
		case 23:
			return integrate_with_boost<adams_bashforth<multistep_order, std::vector<double>>>(tw, t0, tend, x, param);
			break;
		case 24:
			return integrate_with_boost<adams_bashforth_moulton<multistep_order, std::vector<double>>>(tw, t0, tend, x, param);
			break;
		case 25:
			return integrate_with_boost_controlled<runge_kutta_fehlberg78<std::vector<double>>>(tw, t0, tend, x, param);
			break;
		case 26:
			return integrate_with_boost_controlled<runge_kutta_dopri5<std::vector<double>>>(tw, t0, tend, x, param);
			break;
		case 27:
			return integrate_with_boost_controlled<runge_kutta_cash_karp54<std::vector<double>>>(tw, t0, tend, x, param);
			break;
		case 28:
			// hat schon eine Schrittweitensteuerung, deshalb nicht _controlled
			return integrate_with_boost<bulirsch_stoer<std::vector<double>>>(tw, t0, tend, x, param);
			break;
		case 29:
			break;
		default:
			// should never happen
			assert(false);
	}

	// should never happen
	assert(false);
	return -1;
}
#endif


int Butcher::RungeKuttaF(TransWorhpProblem *tw, double &t, double t0, double tend, double *x, double *u1, double *u2, double *param, int startflag) {

	if (!isInit) {
		MyStatus("Butcher", "Butcher-Tableau not initialized! (use explTW::butcherInit())",
		         Status::WARN);
		return -1;
	}

	ph = tw;

	int stepcount=0;

	int ret = r8_rkf45 ( f, n_ode, x, tmp_x.data(), &t, t0, tend, u1, u2, param, &relerr, abserr, startflag, tw->n_ctrl, stepcount );

	if (ret==2 || ret==-2) {
		return stepcount;
	}
	return 0;
}

int Butcher::RungeKutta(TransWorhpProblem *tw, double &t, const double t0, const double tend, double *x, const double *u1, const double *u2, const double *param, double &h) {

	if (!isInit) {
		MyStatus("Butcher", "Butcher-Tableau not initialized! (use explTW::butcherInit())",
		         Status::WARN);
		return -1;
	}

	if (hasStepSize()) {// mit Schrittweitensteuerung

		for (int i=0; i<s; ++i) {

			for (int ode = 0; ode < n_ode; ++ode) {
				if (i>0) {
					double aux = k[0][ode]*a[i][0];
					for (int z = 1;  z < i; ++z) {
						aux += k[z][ode]*a[i][z];
					}
					tmp_x[ode] = x[ode] + aux*h;
				} else {
					tmp_x[ode] = x[ode];
				}
			}

			const double ttt = t+c[i]*h;
			const double pos = (ttt-t0)/(tend-t0);

			for (int l=0; l < n_ctrl; l++) {
				if (tw->solver->twparameter->linInter==1) { // lineare Interpolation der Steuerung
					tmp_u[l] = u1[l] + (u2[l]-u1[l]) * pos;
				} else { // konstante Steuerung
					tmp_u[l] = u1[l];
				}
			}

			tw->ode(&k[i][0], ttt, tmp_x.data(), tmp_u.data(), param);
		}

		for (int ode=0; ode<n_ode; ode++) {
			tmp_ret1[ode] = k[0][ode]*b[0];
			tmp_ret2[ode] = k[0][ode]*b2[0];
			for (int i = 1; i < s; i++) {
				tmp_ret1[ode] += k[i][ode]*b[i];
				tmp_ret2[ode] += k[i][ode]*b2[i];
			}
		}

		double est = 0; // Fehlerabschaetzung
		for (int l = 0; l < n_ode; l++) {
			double tmp = (tmp_ret1[l] - tmp_ret2[l])/h;
			if (tmp < 0) tmp = -tmp;
			if (tmp > est) est = tmp;
		}

		double epsEst;
		if (est == 0) { //dh beide Lsg sind gleich
			epsEst = maxEpsEst;
		} else {
			epsEst = eps/est;
		}

		double factor; // Faktor mit dem die alte SW skaliert wird
		if (epsEst < maxEpsEst && epsEst > minEpsEst) {
			factor = aa*pow(epsEst,1./p);

			if (factor>amax) {
				factor = amax;
			} else if (factor<amin) {
				factor = amin;
			}
		} else {
			if (epsEst >= maxEpsEst) {
				factor = amax;
			} else {
				factor = amin;
			}
		}

		const double hneu = factor * h; // neue Schrittweite

		if (hneu<hmin) { // Schrittweite zu klein
			return -1;

		} else if (est<eps) { // positiver Abbruch
			for (int l=0; l<n_ode; l++) {
				x[l] += tmp_ret1[l]*h;
			}
			t += h;
			h = hneu;

			return 0;

		} else { // Genauigkeit noch nicht erfuellt
			h = hneu;
			return 1;
		}

		//cerr << "Schrittweite: " << h << endl;
	}

	// ohne Schrittweitensteuerung
	else {

		for (int i = 0; i < s; i++) {

			for (int ode = 0; ode < n_ode; ode++) {
				if (i>0) {
					double aux = k[0][ode]*a[i][0];
					for (int z = 1;  z < i; z++) {
						aux += k[z][ode]*a[i][z];
					}
					tmp_x[ode] = x[ode] + aux*h;
				} else {
					tmp_x[ode] = x[ode];
				}
			}

			double ttt = t+c[i]*h;
			double pos = (ttt-t0)/(tend-t0);

			for (int l = 0; l < n_ctrl; l++) {
				if (tw->solver->twparameter->linInter==1) { // lineare Interpolation der Steuerung
					tmp_u[l] = u1[l] + (u2[l]-u1[l]) * pos;
				} else { // konstante Steuerung
					tmp_u[l] = u1[l];
				}
			}

			tw->ode(&k[i][0], ttt, tmp_x.data(), tmp_u.data(), param);
		}

		std::fill(tmp_ret1.begin(),tmp_ret1.end(),0.0);

		for (int i = 0; i < s; i++) {
			for (int l = 0; l < n_ode; l++) {
				tmp_ret1[l] += k[i][l]*b[i];
			}
		}

		for (int l = 0; l < n_ode; l++) {
			x[l] += tmp_ret1[l]*h;
		}

		t += h;

		return 0; // Schritt ist immer "erfolgreich"
	}
}


#if defined(TW_WITH_SUPERLU) || defined(WITH_LAPACK)
int Butcher::ROW(ExplTransWorhp* tw, const double t0, const double tend, double* x, const double* u1, const double* u2, const double* param) {

#ifdef TW_WITH_SUPERLU
	return ROW_SUPERLU(tw,t0,tend,x,u1,u2,param);
#endif

#ifdef WITH_LAPACK
	return ROW_LAPACK(tw,t0,tend,x,u1,u2,param);
#endif
}
#endif

#ifdef TW_WITH_SUPERLU
// nur fuer autonome DGL!!
int Butcher::ROW_SUPERLU(ExplTransWorhp* tw, const double t0, const double tend, double* x, const double* u1, const double* u2, const double* param) {

	// (ein mal) pruefen, ob Ableitungen gesetzt sind
	if (!ROW_DS_check){
		if (!tw->ode_diff2(ROW_DS,0.0,x,tmp_u,param)) {
			cout << "FEHLER: ode_diff2() implementieren!";
			return -1;
		}
	}

	// Schrittweite
	//const double h = 0.01;
	// Anzahl der Schritte
	//const int N = (tend-t0)/(h);

	// Anzahl der Schritte
	const int N = 10;
	// Schrittweite
	const double h = (tend-t0)/N;

	// Anzahl nicht-Null Eintraege der Matrix A bestimmen
	if (nnzOfA == 0) {
		for (int ix = 0; ix < n_ode; ++ix) {
			for (int iy = 0; iy < n_ode; ++iy) {
				if (ROW_DS.check(iy,ix) ||(ix == iy) ) {
					++nnzOfA;
				}
			}
		}
	}
	//cerr << "NNZ: " << nnz << endl;

	//cerr << "Anzahl: " << N << ", Schrittweite: " << h << endl;

	//Speicher holen
	double *rhs1 = new double[n_ode];
	double *rhs2 = new double[n_ode];

	// fuer "ohne struktur"
	//nnz = n_ode*tw->n_ode;

	// Speicher holen fuer A
	int *asub = new int[nnzOfA];
	int *xa = new int[n_ode+1];
	double *aValues = new double[nnzOfA];


	for (int dis = 0; dis < N; ++dis) {

		const double t = t0+h*dis;
		const double pos = (t-t0)/(tend-t0);

		for (int i=0; i<tw->n_ctrl; ++i) {
			if (tw->twparameter->linInter==1) { // lineare Interpolation der Steuerung
				tmp_u[i] = u1[i] + (u2[i]-u1[i]) * pos;
			} else { // konstante Steuerung
				tmp_u[i] = u1[i];
			}
		}

		SuperMatrix A, AC, L, U, B1, B2;
		int *perm_r, *perm_c;
		int *etree;
		int info;
		superlu_options_t options;
		SuperLUStat_t stat;
		GlobalLU_t Glu;

		/* Initialize matrix A. */
		const int n = n_ode;

		tw->ode_diff2(ROW_DS,0.0,x,tmp_u,param);


		// ohne struktur (dense)
		/*
		nnz = n_ode*n_ode;
		for (int ix = 0, i = 0; ix < n_ode; ++ix) {
			xa[ix] = i;
			for (int iy = 0; iy < n_ode; ++iy, ++i) {
				asub[i] = iy;
				a[i] = -yy[0][0]*h*ROW_DS.get(iy,ix);
				if (ix == iy) {
					a[i] += 1.0;
				}
			}
		}
		xa[n_ode] = nnz;
		*/

		// mit struktur
		for (int ix = 0, i = 0; ix < n_ode; ++ix) {
			xa[ix] = i;
			for (int iy = 0; iy < n_ode; ++iy) {
				if (ROW_DS.check(iy,ix)){
					asub[i] = iy;
					aValues[i] = -h*gamma[0][0]*ROW_DS.get(iy,ix);
					if (ix == iy) {
						aValues[i] += 1.0;
					}
					++i;
				}
				else if (ix == iy) {
					asub[i] = iy;
					aValues[i] = 1.0;
					++i;
				}
			}
		}
		xa[n] = nnzOfA;


		// A Matrix erstellen -> (1 - h \gamma T) (bleibt immer gleich)
		dCreate_CompCol_Matrix(&A, n, n, nnzOfA, aValues, asub, xa, SLU_NC, SLU_D, SLU_GE);

		// print matrix
		//dPrint_CompCol_Matrix("A", &A);

		// Dynamik auswerten
		tw->ode(tmp_x.data(),t0+h*c[0],x,tmp_u,param);

		// Werte rechte Seite (fuer k1)
		for (int i = 0; i < n_ode; ++i) {
			rhs1[i] = tmp_x[i];
		}

		// rechte Seite erstellen (dense) (fuer k1)
		dCreate_Dense_Matrix(&B1, n, 1, rhs1, n, SLU_DN, SLU_D, SLU_GE);

		// print matrix
		//dPrint_Dense_Matrix("B",&B);

		// Speicher holen
		if ( !(perm_r = intMalloc(n)) ) ABORT("Malloc fails for perm_r[].");
		if ( !(perm_c = intMalloc(n)) ) ABORT("Malloc fails for perm_c[].");

		// Set the default input options.
		set_default_options(&options);
		options.ColPerm = NATURAL;

		// Initialize the statistics variables.
		StatInit(&stat);

		// loese LGS fuer k1
		//dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);

		/////////////////////////
		// erst faktorisieren

		const int panel_size = sp_ienv(1);
		const int relax = sp_ienv(2);

		if ( !(etree = intMalloc(n)) ) ABORT("Malloc fails for etree[].");

		get_perm_c(options.ColPerm, &A, perm_c);

		sp_preorder(&options, &A, perm_c, etree, &AC);

		// faktorisieren
		dgstrf(&options, &AC, relax, panel_size, etree, NULL, 0, perm_c, perm_r, &L, &U, &Glu, &stat, &info);


		// loese LGS fuer k1 (nutze Faktorisierung)
		dgstrs(NOTRANS,&L,&U,perm_c,perm_r, &B1,&stat,&info);

		////////////////////////

		// Loesung ausgeben
		//dPrint_Dense_Matrix("loesung ",&B);

		// speichere k1
		double *k1 = (double*) ((DNformat*) B1.Store)->nzval;

		// x extraploieren
		for (int i = 0; i < n_ode; ++i) {
			tmp_ret1[i] = x[i] + h*a[1][0]*k1[i];
		}

		// Dynamik auswerten (eig autonom)
		tw->ode(tmp_x.data(),t0+h*c[1],tmp_ret1,tmp_u,param);

		// Werte rechte Seite (fuer k2)
		for (int i = 0; i < n_ode; ++i) {
			double sum = 0;
			for (int j = 0; j < n_ode; ++j) {
				sum += ROW_DS.get(i,j)*k1[j];
			}
			rhs2[i] = tmp_x[i] + h*gamma[1][0]* sum;
		}

		// rechte Seite erstellen (dense) (fuer k2)
		dCreate_Dense_Matrix(&B2, n, 1, rhs2, n, SLU_DN, SLU_D, SLU_GE);

		// print matrix
		//dPrint_Dense_Matrix("B2",&B2);

		// loese LGS fuer k2 (nutze Faktorisierung)
		dgstrs(NOTRANS,&L,&U,perm_c,perm_r, &B2,&stat,&info);

		//dPrint_Dense_Matrix("loesung ",&B2);

		// speichere k2
		double *k2 = (double*) ((DNformat*) B2.Store)->nzval;

		for (int i = 0; i < n_ode; ++i) {
			x[i] = x[i] + h*(b[0]*k1[i]+b[1]*k2[i]);
		}
		//cerr << "approx: " << x[0] << " " << x[1] << endl;

		// Speicher freigeben
		SUPERLU_FREE (perm_r);
		SUPERLU_FREE (perm_c);

		SUPERLU_FREE (etree);

		//Destroy_CompCol_Matrix(&A);
		Destroy_SuperMatrix_Store(&A);
		Destroy_CompCol_Permuted(&AC);

		Destroy_SuperMatrix_Store(&B1);
		Destroy_SuperMatrix_Store(&B2);
		Destroy_SuperNode_Matrix(&L);
		Destroy_CompCol_Matrix(&U);
		StatFree(&stat);
	}

	// Speicher freigeben
	delete[] rhs1;
	delete[] rhs2;

	delete[] asub;
	delete[] xa;
	delete[] aValues;

	return N;
}
#endif

#ifdef WITH_LAPACK
int Butcher::ROW_LAPACK(ExplTransWorhp* tw, const double t0, const double tend, double* x, const double* u1, const double* u2, const double* param) {

	// (ein mal) pruefen, ob Ableitungen gesetzt sind
	if (!ROW_DS_check){
		if (!tw->ode_diff2(ROW_DS,0.0,x,tmp_u,param)) {
			cout << "FEHLER: ode_diff2() implementieren!";
			return -1;
		}
	}

	// Schrittweite
	//const double h = 0.01;
	// Anzahl der Schritte
	//const int N = (tend-t0)/(h);

	// Anzahl der Schritte
	const int N = 10;
	// Schrittweite
	const double h = (tend-t0)/N;

	int n = n_ode;
	int nrhs = 1;

	// LAPACK kann nur dense
	nnzOfA = n*n;

	//cerr << "Anzahl: " << N << ", Schrittweite: " << h << endl;

	//Speicher holen
	double *rhs1 = new double[n];
	double *rhs2 = new double[n];


	// Speicher holen fuer A
	double *A = new double[nnzOfA];

	int *ipiv = new int[n];


	for (int dis = 0; dis < N; ++dis) {

		const double t = t0+h*dis;
		const double pos = (t-t0)/(tend-t0);

		for (int i=0; i<tw->n_ctrl; ++i) {
			if (tw->twparameter->linInter==1) { // lineare Interpolation der Steuerung
				tmp_u[i] = u1[i] + (u2[i]-u1[i]) * pos;
			} else { // konstante Steuerung
				tmp_u[i] = u1[i];
			}
		}

		// Locals
		int info;

		int lda = n;
		int ldb = n;

		tw->ode_diff2(ROW_DS,0.0,x,tmp_u,param);


		// ohne struktur (dense)
		for (int ix = 0, i = 0; ix < n; ++ix) {
			for (int iy = 0; iy < n; ++iy, ++i) {

				A[i] = -gamma[0][0]*h*ROW_DS.get(iy,ix);
				if (ix == iy) {
					A[i] += 1.0;
				}
			}
		}


		// Dynamik auswerten
		tw->ode(tmp_x.data(),t0+h*c[0],x,tmp_u,param);

		// Werte rechte Seite (fuer k1)
		for (int i = 0; i < n; ++i) {
			rhs1[i] = tmp_x[i];
		}

		// loese LGS fuer k1
		//dgesv_(n, nrhs, A, lda, ipiv, rhs1, ldb, info);

		// loese LGS fuer k1 (mit Faktorisierung)
		dgetrf_(n, n, A, lda, ipiv, info);
		dgetrs_("No transpose", n, nrhs, A, lda, ipiv, rhs1, ldb, info );


		// speichere k1
		double *k1 = rhs1;

		// x extraploieren
		for (int i = 0; i < n; ++i) {
			tmp_ret1[i] = x[i] + h*a[1][0]*k1[i];
		}

		// Dynamik auswerten (eig autonom)
		tw->ode(tmp_x.data(),t0+h*c[1],tmp_ret1,tmp_u,param);

		// Werte rechte Seite (fuer k2)
		for (int i = 0; i < n; ++i) {
			double sum = 0.0;
			for (int j = 0; j < n_ode; ++j) {
				sum += ROW_DS.get(i,j)*k1[j];
			}
			rhs2[i] = tmp_x[i] + h*gamma[1][0]* sum;
		}

		// loese LGS fuer k2
		//dgesv_(n, nrhs, A, lda, ipiv, rhs2, ldb, info);

		// loese LGS fuer k2 (mit Faktorisierung)
		dgetrs_("No transpose", n, nrhs, A, lda, ipiv, rhs2, ldb, info );

		// speichere k2
		double *k2 = rhs2;

		for (int i = 0; i < n; ++i) {
			x[i] = x[i] + h*(b[0]*k1[i]+b[1]*k2[i]);
		}
		//cerr << "approx: " << x[0] << " " << x[1] << endl;

	}

	// Speicher freigeben
	delete[] rhs1;
	delete[] rhs2;

	delete[] A;
	delete[] ipiv;

	return N;
}
#endif

}
