/*----------------------------------------------------------------
 *
 *  Example: Betts Testset
 *
 *----------------------------------------------------------------*/

#include "problems.h"

#include "TWfolder.h"
#include "TWproblem.h"
#include "TransWORHP.h"

#include <memory>
#include <chrono>

extern "C"
{
	void setiprob_(int &n);
}

using ODE = void (*)(int&,double&,const double*,int&,const double*,const int&,double*,const int&,int&);
using INI = void (*)(int&,int&,int&,int&,int*,int*,int&,int&,int&,int&,int&,const int&,double*,double*,double*,char*,const int&,double*,double*,double*,double*,double*,char*,const int&,double*,double*,char*,const int&,double*,int*,char*,int&);

struct Problem {
	
	std::string probName;
	ODE dynamic;
	INI initalisation;
	int probNo;
	
	
	Problem(int _probNo, std::string name, ODE _dynamic, INI _initalisation) :
		probName(std::move(name)), dynamic(_dynamic), initalisation(_initalisation), probNo(_probNo) {
	}
};

struct TestCollection {

private:
	std::vector<Problem> tests;

public:
	TestCollection() {
		tests.reserve(51);
		tests.emplace_back(178, "aqua01", aquade_, aquain_); //aqua01 //  0
		tests.emplace_back(151, "arao01", araode_, araoin_); //arao01 //  1
		tests.emplace_back(152, "arao02", araode_, araoin_); //arao02 //  2
		tests.emplace_back(222, "medi02", medide_, mediin_); //medi02 //  3
		tests.emplace_back(224, "medi04", medide_, mediin_); //medi04 //  4
		tests.emplace_back(226, "medi06", medide_, mediin_); //medi06 //  5
		tests.emplace_back(176, "pndl02", pndlde_, pndlin_); //pndl02 //  6
		tests.emplace_back( 61, "brgr01", brgrde_, brgrin_); //brgr01 //  7
		tests.emplace_back(127, "ashr01", ashrde_, ashrin_); //ashr01 //  8
		tests.emplace_back(143, "ashr02", ashrde_, ashrin_); //ashr02 //  9
		tests.emplace_back(144, "ashr03", ashrde_, ashrin_); //ashr03 // 10
		tests.emplace_back(145, "ashr04", ashrde_, ashrin_); //ashr04 // 11
		tests.emplace_back(153, "ashr05", ashrde_, ashrin_); //ashr05 // 12
		tests.emplace_back(274, "ashr06", ashrde_, ashrin_); //ashr06 // 13
		tests.emplace_back(105, "brac01", bracde_, bracin_); //brac01 // 14
		tests.emplace_back(131, "chmr01", chmrde_, chmrin_); //chmr01 // 15
		tests.emplace_back(132, "chmr02", chmrde_, chmrin_); //chmr02 // 16
		tests.emplace_back(133, "chmr03", chmrde_, chmrin_); //chmr03 // 17
		tests.emplace_back(134, "chmr04", chmrde_, chmrin_); //chmr04 // 18
		tests.emplace_back(135, "chmr05", chmrde_, chmrin_); //chmr05 // 19
		tests.emplace_back(136, "chmr06", chmrde_, chmrin_); //chmr06 // 20
		tests.emplace_back(137, "chmr07", chmrde_, chmrin_); //chmr07 // 21
		tests.emplace_back(138, "chmr08", chmrde_, chmrin_); //chmr08 // 22
		tests.emplace_back(139, "chmr09", chmrde_, chmrin_); //chmr09 // 23
		tests.emplace_back(140, "chmr10", chmrde_, chmrin_); //chmr10 // 24
		tests.emplace_back(183, "cran01", crande_, cranin_); //cran01 // 25
		tests.emplace_back(190, "ffrb01", ffrbde_, ffrbin_); //ffrb01 // 26
		tests.emplace_back(123, "hang01", hangde_, hangin_); //hang01 // 27
		tests.emplace_back(124, "hang02", hangde_, hangin_); //hang02 // 28
		tests.emplace_back(273, "hang03", hangde_, hangin_); //hang03 // 29
		tests.emplace_back(191, "jshi01", jshide_, jshiin_); //jshi01 // 30
		tests.emplace_back(192, "jshi02", jshide_, jshiin_); //jshi02 // 31
		tests.emplace_back(193, "lnht01", lnhtde_, lnhtin_); //lnht01 // 32
		tests.emplace_back(194, "lnht02", lnhtde_, lnhtin_); //lnht02 // 33
		tests.emplace_back(  7, "lnts05", lntsde_, lntsin_); //lnts05 // 34
		tests.emplace_back(271, "lnts13", lntsde_, lntsin_); //lnts13 // 35
		tests.emplace_back(130, "lowt01", lowtde_, lowtin_); //lowt01 // 36
		tests.emplace_back(294, "rbrm01", rbrmde_, rbrmin_); //rbrm01 // 37
		tests.emplace_back(211, "tb2s01", tb2sde_, tb2sin_); //tb2s01 // 38
		tests.emplace_back(180, "tran01", trande_, tranin_); //tran01 // 39
		tests.emplace_back(180, "zrml01", zrmlde_, zrmlin_); //zrml01 // 40
		tests.emplace_back(164, "alpr01", alprde_, alprin_); //alpr01 // 41
		tests.emplace_back(129, "brac01", bracde_, bracin_); //brac01 // 42
		tests.emplace_back(148, "heat01", heatde_, heatin_); //heat01 // 43
		//tests.emplace_back(150,"heat02", heatde_, heatin_); //heat02 // 44
		//tests.emplace_back(128,"kplr01", kplrde_, kplrin_); //kplr01 // 45
		tests.emplace_back(179, "pndl02", pndlde_, pndlin_); //pndl02 // 46
		tests.emplace_back(  1, "qlin01", qlinde_, qlinin_); //qlin01 // 47
		tests.emplace_back(  2, "qlin02", qlinde_, qlinin_); //qlin02 // 48
		tests.emplace_back(146, "qlin03", qlinde_, qlinin_); //qlin03 // 49
		//tests.emplace_back(147, "qlin04", qlinde_, qlinin_); //qlin04 // 50
	}
	
	size_t size() const {
		return tests.size();
	}
	
	const Problem& getProblem(int no) const {
		return tests[no];
	}
};


class TestSet : public tw::TransWorhpProblem {

private:
	
	int phase;
	int maxOrMin;
	int initType;
	
	ODE dynamic;
	
	double *XLB;
	double *XUB;
	
	double *PLB;
	double *PUB;
	
	double *CLB;
	double *CUB;
	
	double *X0;
	double *X1;
	
	double *P0;
	
	int *ITERM;
	double *COEF;
	
	char *label;
	
public:
	
	TestSet(char *title, char *_label, int dis, int ode, int ctrl, int param, int neben, int integral,
		int _phase, int MAXMIN,
		double *YLB, double *YUB,
		double *_PLB, double *_PUB,
		double *_CLB, double *_CUB,
		double *Y0, double *Y1,
		int ini,
		double *_P0,
		int *_ITERM, double *_COEF,
		ODE dyn) :
		TransWorhpProblem(tw::TWdimension(title,dis,ode,ctrl,param,0,neben,integral)),
		phase(_phase),
		maxOrMin(MAXMIN),
		initType(ini),
		dynamic(dyn),
		XLB(YLB), XUB(YUB),
		PLB(_PLB), PUB(_PUB),
		CLB(_CLB), CUB(_CUB),
		X0(Y0), X1(Y1),
		P0(_P0),
		ITERM(_ITERM), COEF(_COEF),
		label(_label) {
		
		// Startschaetzung-Typ abfragen
		if (initType > 2 || initType == 0) {
			std::cout << "Fehler initType " << initType << std::endl;
			if (initType != 6) {
				exit(1);
			}
		}
	}
	
	void localinit() override {
		if (XLB[0] == XUB[0] && XLB[2] == XUB[2]) {
			solver->LinearTimeAxis(XLB[0],XLB[2]);
		} else if (n_param > 0) {
			std::cout << "freetime" << std::endl;
			//cout << XLB[0] << " " << XUB[0] << " " << XLB[2] << " " << XUB[2] << endl;
			freetime = true;
		} else {
			std::cout << XLB[0] << " " << XUB[0] << " " << XLB[2] << " " << XUB[2] << std::endl;
			std::cout << "Fehler" << std::endl;
			exit(1);
		}
	}
	
	
	std::string GetXTitle(int d) override {
		int counter = 0;
		const int index = 80+(d+1)*79;
		for (int i = 0; i < 80; i++) {
			if (label[index-i] == ' ' || label[index-i] == 0) {
				counter++;
			} else {
				break;
			}
		}
		return std::string(&label[80+d*80],80-counter);
	}

	std::string GetUTitle(int d) override {
		int counter = 0;
		const int index = 80+n_ode*80+(d+1)*79;
		for (int i = 0; i < 80; i++) {
			if (label[index-i] == ' ' || label[index-i] == 0) {
				counter++;
			} else {
				break;
			}
		}
		return std::string(&label[80+n_ode*80+d*80],80-counter);
	}
	
	double obj() override {
		
		double ret = 0;
		
		for (int i = 0; i < 60; i+=4) {
			if (ITERM[i] == 0) { // ist Zielfunktion
				if (ITERM[i+1] == phase) { // richtige Phase
					if(ITERM[i+2] == +1) { // Auswertung am Ende
						if (ITERM[i+3] != 0) {
							ret += (-1)*maxOrMin * x(n_dis-1,abs(ITERM[i+3])-1);
						} else {
							if (freetime) {
								ret += p(0);
							} else {
								std::cout << "Fehler obj" << std::endl;
								exit(1);
							}
						}
					}
				}
			} else if (ITERM[i] == -1) {
				break;
			}
		}

		return ret;
	}

	
	bool obj_structure(tw::DiffStructure &s) override {
	  
		bool set = false;
	  
		for (int i = 0; i < 60; i+=4) {
			if (ITERM[i] == 0) { // ist Zielfunktion
				if (ITERM[i+1] == phase) { // richtige Phase
					if(ITERM[i+2] == +1) { // Auswertung am Ende
						if (ITERM[i+3] != 0) {
							s(0, x_index(n_dis-1,abs(ITERM[i+3])-1));
							set = true;
						} else {
							if (freetime) {
								s(0, p_index(0));
								set = true;
							} else {
								std::cout << "Fehler obj" << std::endl;
								exit(1);
							}
						}
					}
				}
			} else if (ITERM[i] == -1) {
				break;
			}
		}
		
		return set;
	}
	
	
	void integral(double *f, double t, const double *x, const double *u, const double *p) override {
		
		int N = n_ode+n_ctrl;
		std::vector<double> aux;
		aux.reserve(N);
		std::copy(x,x+n_ode,aux.begin());
		std::copy(u,u+n_ctrl,aux.begin()+n_ode);
		std::vector<double> aux2(n_ode+n_neben+n_integral,0);
		int iferr = 0;
		
		dynamic(phase,t,aux.data(),N,p,n_param,aux2.data(),n_ode+n_neben+n_integral,iferr);
		
		for (int i = 0; i < n_integral; i++) {
			f[i] = (-1)*maxOrMin * COEF[n_neben+i] * aux2.at(n_ode+n_neben+i);
		}
	}

	void ode(double *dx, double t, const double *x, const double *u, const double *p) override {
		
		int N = n_ode+n_ctrl;
		std::vector<double> aux;
		aux.reserve(N);
		std::copy(x,x+n_ode,aux.begin());
		std::copy(u,u+n_ctrl,aux.begin()+n_ode);
		std::vector<double> aux2(n_ode+n_neben+n_integral,0);
		int iferr = 0;
		
		if (!freetime) {
			dynamic(phase,t,aux.data(),N,p,n_param,aux2.data(),n_ode+n_neben+n_integral,iferr);
			for (int i = 0; i < n_ode; ++i) {
				dx[i] = aux2[i];
			}
		} else {
			dynamic(phase,t,aux.data(),N,&p[1],n_param,aux2.data(),n_ode+n_neben+n_integral,iferr);
			for (int i = 0; i < n_ode; ++i) {
				dx[i] = aux2[i] * p[0];
			}
		}
	}
	
	
	void neben(double *c, double t, const double *x, const double *u, const double *p) override {
		
		int N = n_ode+n_ctrl;
		std::vector<double> aux;
		aux.reserve(N);
		std::copy(x,x+n_ode,aux.begin());
		std::copy(u,u+n_ctrl,aux.begin()+n_ode);
		std::vector<double> aux2(n_ode+n_neben+n_integral,0);
		int iferr = 0;
		
		dynamic(phase,t,aux.data(),N,&p[1],n_param,aux2.data(),n_ode+n_integral,iferr);
		
		for (int i = 0; i < n_neben; ++i) {
			c[i] = COEF[i] * aux2.at(n_ode+i);
		}
	}
	
	void neben_boundary(double *c_low, double *c_upp) {
		
		for (int i = 0; i < n_neben; ++i) {
			c_low[i] = CLB[i];
			c_upp[i] = CUB[i];
		}
	}
	
	void x_init(double *x, int i, int dis) override {
		if (i == 0) {
			for (int k = 0; k < n_ode; k++) {
				x[k] = X0[1+k];
			}
		} else if (i == dis-1) {
			for (int k = 0; k < n_ode; k++) {
				x[k] = X1[1+k];
			}
		} else {
			for (int k = 0; k < n_ode; k++) {
				x[k] = i * (X1[1+k]-X0[1+k])/(dis-1) + X0[1+k];
			}
		}
	}
	
	void u_init(double *u, int i, int dis) override {
		const int index = 1 + n_ode;
		if (i == 0) {
			for (int k = 0; k < n_ctrl; k++) {
				u[k] = X0[index + k];
			}
		} else if (i == dis-1) {
			for (int k = 0; k < n_ctrl; k++) {
				u[k] = X1[index + k];
			}
		} else {
			for (int k = 0; k < n_ctrl; k++) {
				u[k] = i * (X1[index + k]-X0[index + k])/(dis-1) + X0[index + k];
			}
		}
	}
	
	void p_init(double *p) override {
		if (!freetime) {
			for (int i = 0; i < n_param; i++) {
				p[i] = P0[i];
			}
		} else {
			p[0] = X1[0];
			for (int i = 1; i < n_param; i++) {
				p[i] = P0[i-1];
			}
		}
	}
	
	
	void x_boundary(double *x_low, double *x_upp) override {
		
		for (int i = 0; i < n_ode; i++) {
			x_low[i] = XLB[3 + 1 + i*3];
			x_upp[i] = XUB[3 + 1 + i*3];
		}
	}
	
	void u_boundary(double *u_low, double *u_upp) override {
		
		for (int i = 0; i < n_ctrl; ++i) {
			u_low[i] = XLB[3 + n_ode*3 + 1 + i*3];
			u_upp[i] = XUB[3 + n_ode*3 + 1 + i*3];
		}
	}
	
	void p_boundary(double *p_low, double *p_upp) override {
		
		if (!freetime) {
			for (int i = 0; i < n_param; i++) {
				p_low[i] = PLB[i];
				p_upp[i] = PUB[i];
			}
		} else {
			p_low[0] = XLB[2];
			p_upp[0] = XUB[2];
			for (int i = 1; i < n_param; i++) {
				p_low[i] = PLB[i-1];
				p_upp[i] = PUB[i-1];
			}
		}
	}

	void var_boundary(double *x_low, double *x_upp) override {
		
		for (int i = 0; i < n_ode; i++) {
			x_low[x_index(0,i)] = XLB[3 + 0 + i*3];
			x_upp[x_index(0,i)] = XUB[3 + 0 + i*3];
			
			x_low[x_index(n_dis-1,i)] = XLB[3 + 2 + i*3];
			x_upp[x_index(n_dis-1,i)] = XUB[3 + 2 + i*3];
		}
		
		for (int i = 0; i < n_ctrl; i++) {
			x_low[u_index(0,i)] = XLB[3 + n_ode*3 + i*3];
			x_upp[u_index(0,i)] = XUB[3 + n_ode*3 + i*3];
			
			x_low[u_index(n_dis-1,i)] = XLB[3 + n_ode*3 + 2 + i*3];
			x_upp[u_index(n_dis-1,i)] = XUB[3 + n_ode*3 + 2 + i*3];
		}
	}
	
	/*
	bool step() {
		static int counter = 0;
		counter++;
		
		if (counter == 2)
			return false;
		
		return true;
	}
	*/
};

/////////////////////////////////////////////////////////////////////////////

struct bettsOut {
	
	std::string probName;

	/** #NLP-Schritte */
	int NLPsteps;
	/** #QP-Schritte */
	int QPsteps;

	/** Optimalitaet */
	double opt;
	/** Zulaessigkeit */
	double con;
	/** Wert Zielfunktion*/
	double obj;
	/** Rechenzeit */
	double time;

	/** erfolgreich? */
	bool isOptimal;
	
	friend std::ostream& operator<<(std::ostream& os, const bettsOut& b);
};

std::ostream& operator<<(std::ostream& os, const bettsOut& b) {
	if (b.isOptimal) {
		os << "\033[32m";
	} else {
		os << "\033[31m";
	}
	os << b.probName << "\t" << b.isOptimal << "\t" << std::setw(4) << b.NLPsteps << "\t" << std::setw(4) << b.QPsteps << "\t" << std::setw(12) << b.opt << "\t" << std::setw(12) << b.con << "\t" << std::setw(12) << b.obj << "\t" << b.time;
	os << "\033[30m";
	return os;
}

int main(int argv, char* argc[]) {


TestCollection myTestCollection;

std::vector<bettsOut> output;
output.reserve(myTestCollection.size());

std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();

bool debug = false;

for (size_t k = 0; k < myTestCollection.size(); k++)
{
	std::cout << "####################### " << std::setw(2) << k << "#######################" << std::endl;

	int IPHASE = 1;			//The current phase number, k.
	int NPHS = 1;			//Number of phases in the problem. (default = 1)
	int METHOD;			//Discretization/integration method, used in conjunction with NSTG.
	int NSTG;			//Discretization stages
	std::vector<int> NCF(6);		//Number of continuous functions
	std::vector<int> NPF(2);		//Number of computed point functions.
	int NPV = 0;			//Number of parameters, p.
	int NAV = 0;			//Number of algebraic or control variables, u.
	int NGRID;			//Number of grid points in initial grid.
	int INIT = 1;			//Initial guess type for internal states
	int MAXMIN = -1;		//Optimization Type Flag.
	const int MXPARM = 5;		//Maximum number of parameters per phase.
	double P0[MXPARM];		//Array of length NPV containing the initial guess for parameters, p.
	double PLB[MXPARM];		//Array of length NPV containing the lower bounds for parameters, p.
	double PUB[MXPARM];		//Array of length NPV containing the upper bounds for parameters, p.
	char plbl[MXPARM*80];		//Symbol and description for parameters, p.
	const int MXSTAT = 60;		//Maximum number of dynamic variables per phase.
	double Y0[MXSTAT+1];		//Real array of length NCF(1)+NAV+1 containing an initial guess for the dynamic variables at the beginning of the phase.
	double Y1[MXSTAT+1];		//Real array of length NCF(1)+NAV+1 containing an initial guess for the dynamic variables at the end of the phase.
	double YLB[(MXSTAT+1)*3];	//Real array of length NCF(1)+NAV+1 containing the lower bounds for the independent and dynamic variables in the interior and boundaries of the phase.
	double YUB[(MXSTAT+1)*3];	//Real array of length NCF(1)+NAV+1 containing the upper bounds for the independent and dynamic variables in the interior and boundaries of the phase.
	double STSKL[(MXSTAT+MXPARM+1)*2];//Real array of length 2*(NCF(1)+NAV+NPV+1) containing scale information
	std::vector<char> stlbl(5120,' ');	//Character array of length NCF(1)+NAV+1 containing symbol and description for independent variable t and dynamic variables, z.
	const int MXPCON = 20;		//Maximum number of user defined constraints per phase.
	double CLB[MXPCON];		//Array containing the lower bounds for constraints.
	double CUB[MXPCON];		//Array containing the upper bounds for constraints.
	char clbl[1024];		//Array containing the symbol and description for constraints.
	const int MXTERM = 50;		//Maximum number of terms allowed in ITERM per phase.
	std::vector<double> COEF(MXTERM);	//Multipliers for user defined terms.
	std::vector<int> ITERM(MXTERM*4,-1);	//Integer array containing a description of the user defined term j used in forming general constraints and/or objective.
	char TITLE[3*60];		//Array containing the titles describing the problem.
	int IER;			//Input Evaluation Error Flag.
	
	
	std::fill(std::begin(YLB),std::end(YLB),-1e200);
	std::fill(std::begin(YUB),std::end(YUB),+1e200);
	
	std::fill(std::begin(CLB),std::end(CLB),-1e200);
	std::fill(std::begin(CUB),std::end(CUB),+1e200);
	
	tw::TWparameter twparameter("transworhp.xml");
	std::map<std::string,std::string> args = twparameter.Arguments(argv,argc);
	tw::TWfolder folder(&twparameter,0);
	
	//int k = ToInt(args["no"]);
	
	Problem myTest = myTestCollection.getProblem(k);
	
	// Problem Nummer im Fortran common-Block setzen
	setiprob_(myTest.probNo);
	
	myTest.initalisation(IPHASE,NPHS,METHOD,NSTG,NCF.data(),NPF.data(),NPV,NAV,NGRID,
		INIT,MAXMIN,MXPARM,P0,PLB,PUB,plbl,
		MXSTAT,Y0,Y1,YLB,YUB,STSKL,stlbl.data(),MXPCON,CLB,CUB,
		clbl,MXTERM,COEF.data(),ITERM.data(),TITLE,IER);
	
	if (debug) {
		std::cout << "(" << YLB[0] << "," << YUB[0] << ") - " << "(" << YLB[2] << "," << YUB[2] << ")" << std::endl;
		for (int i = 3; i < 3+(NCF[0]+NAV)*3; i+=3) {
			std::cout << "(" << YLB[i] << "," << YUB[i] << ") - ";
			std::cout << "(" << YLB[i+1] << "," << YUB[i+1] << ") - ";
			std::cout << "(" << YLB[i+2] << "," << YUB[i+2] << ")" << std::endl;
		}
		//exit(1);
	}
	
	
	// Endzeit ist frei?
	if (YLB[0] == YUB[0] && YLB[2] < YUB[2]) {
		NPV++; // #param erhoehen
	}
	
	// Titel-Text beschneiden
	int counter = 0;
	for (int i = 0; i < 60; i++) {
		if (TITLE[59-i] == ' ' || TITLE[59-i] == 0) {
			counter++;
		} else {
			break;
		}
	}
	TITLE[60-counter] = 0;
	
	TestSet ph(TITLE,stlbl.data(),twparameter.NDIS,NCF[0],NAV,NPV,NCF[1],NCF[2],
		IPHASE,MAXMIN,
		YLB,YUB,
		PLB,PUB,
		CLB,CUB,
		Y0,Y1,
		INIT,
		P0,
		ITERM.data(),COEF.data(),
		myTest.dynamic);
	
	ph.setSolver(&twparameter);
	folder.Add(&ph);
	
	std::unique_ptr<tw::Viewer> viewer;
	if (twparameter.PLOT) viewer = std::unique_ptr<tw::Viewer>(new tw::Viewer(&twparameter));
	
	folder.Init();
	folder.Init(viewer.get());
	
	if (MAXMIN == 0) { //Zulaessigkeit
		folder.worhp_p.FeasibleOnly = true;
	}
	
	std::chrono::steady_clock::time_point t_start = std::chrono::steady_clock::now();
	folder.Loop();
	std::chrono::steady_clock::time_point t_end = std::chrono::steady_clock::now();
	
	bettsOut out;
	
	if (MAXMIN == 0) {
		out.isOptimal = (folder.worhp_c.status==FeasibleSolution);
	} else {
		out.isOptimal = (folder.worhp_c.status==OptimalSolution);
	}
	out.probName = myTest.probName;
	out.NLPsteps = folder.worhp_w.MajorIter;
	out.QPsteps = folder.worhp_w.MinorIterOverall;
	out.opt = folder.worhp_w.NormMax_DL;
	out.con = folder.worhp_w.NormMax_CV;
	out.obj = folder.worhp_o.F/folder.worhp_w.ScaleObj;
	out.time = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
	
	output.push_back(std::move(out));
	
	if (debug) {
		std::cout << "IPHASE         " << IPHASE << std::endl;
		std::cout << "NPHS           " << NPHS << std::endl;
		std::cout << "METHOD         " << METHOD << std::endl;
		std::cout << "NSTG           " << NSTG << std::endl;
		std::cout << "NCF[0] (ODE)   " << NCF[0] << std::endl;
		std::cout << "NCF[1] (DAE)   " << NCF[1] << std::endl;
		std::cout << "NCF[2] (QUAD)  " << NCF[2] << std::endl;
		std::cout << "NCF[3] (AUX)   " << NCF[3] << std::endl;
		std::cout << "NCF[4] (DISC)  " << NCF[4] << std::endl;
		std::cout << "NCF[5] (DELAY) " << NCF[5] << std::endl;
		std::cout << "NPF[0] (POINT) " << NPF[0] << std::endl;
		std::cout << "NPF[1] (POINT) " << NPF[1] << std::endl;
		std::cout << "NPV            " << NPV << std::endl;
		std::cout << "NAV            " << NAV << std::endl;
		std::cout << "NGRID          " << NGRID << std::endl;
		std::cout << "INIT           " << INIT << std::endl;
		std::cout << "MAXMIN         " << MAXMIN << std::endl;
		std::cout << "MXPARM         " << MXPARM << std::endl;
		
		std::cout << std::endl;
		for (int i = 0; i < MXTERM*4; i++) {
			std::cout << "ITERM[" << i << "] " << ITERM[i] <<  std::endl;
		}
		std::cout << "......" << std::endl;
		for (int i = 0; i < MXTERM; i++) {
			std::cout << "COEF " << COEF[i] << std::endl;
		}
		
		std::cout << "...Startschaetzung..." << std::endl;
		for (int i = 0; i < NCF[0]+NAV+1; i++) {
			std::cout << Y0[i] << "\t" << Y1[i] << std::endl;
		}
		std::cout << "...Startschaetzung.p." << std::endl;
		for (int i = 0; i < NPV; i++) {
			std::cout << P0[i] << std::endl;
		}
		
		std::cout << "...CON Schranken..." << std::endl;
		for (int i = 0; i < NCF[1]; i++) {
			std::cout << CLB[i] << " " << CUB[i] << std::endl;
		}
	}
}
	
	std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
	
	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
	std::cout << "It took me " << time_span.count() << " seconds.";
	
	std::cout << std::endl;
	std::cout << "OUTPUT " << output.size() << std::endl;
	std::cout << "\t" << "name" << "\t" << "opti" << "\t" << std::setw(4) << "SQP" << "\t" << std::setw(4) << "QP" << "\t" << std::setw(12) << "kkt" << "\t" << std::setw(12) << "con" << "\t" << std::setw(12) << "obj" << "\t" << "time" << std::endl;
	for (size_t i = 0; i < output.size(); i++) {
		std::cout << i << ":\t" << output[i] << std::endl;
	}
	
	int ave = 0;
	for (auto ele : output) {
		if(ele.isOptimal) ave++;
	}
	std::cout << "I solved " << ave << " of " << output.size() << " problems." << std::endl;
	
	
	return 0;
}

