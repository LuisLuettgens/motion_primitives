/*----------------------------------------------------------------
 *
 *  Example: optimales Fischen
 *
 *----------------------------------------------------------------*/

#ifdef WIN32
#include "windows.h"
#endif

#include "TransWORHP.h"

using namespace std;

class FischenPhase : public TransWorhp {
public:

	double r;
	double k;
	double delta;
	double c;
	double pp;

	FischenPhase(int dis) : TransWorhp("optimales Fischen", dis,2,1,0,0,0) {
	
	
		r = 11;
		k = 5;
		delta = 0.1;
		c = 1;
		pp = 2;
	}

	string GetXTitle(int d) override {
		if (d==0) return "x_1";
		if (d==1) return "aux";
	}

	string GetUTitle(int d) override {
		if (d==0) return "u";
	}
	
	void x_init(double *x, int i, int dis) override {
		x[0] = 4;
		x[1] = 0;
	}
	
	void p_init(double *p) override {
		
	}

	void u_init(double *u, int i, int dis) override {
		u[0] = 0.3;
	}
	
	double obj() override {
		return -x(n_dis-1,1);
	}

	bool obj_structure(DiffStructure &s) override {
	
		s(0,x_index(n_dis-1,1));
		return true;
	}
	
	bool obj_diff(DiffStructure &s) override {
	
		s(0,x_index(n_dis-1,1)) = -1;
		return true;
	}
	

	void ode(double *dx, double t, const double *x, const double *u, const double *p) override {
		dx[0] = (r*x[0]*(1.-x[0]/k) - u[0]);
		dx[1] = (exp(-delta*t)*(pp-c/x[0])*u[0]);
	}

	bool ode_structure(DiffStructure &s) override {
		return false;
	}

	bool ode_diff(DiffStructure& s, double t, const double* x, const double* u, const double* p) override {
		return false;
	}
	
	void x_boundary(double *x_low, double *x_upp) override {
		x_low[0] = 0;
		x_low[1] = 0;
	}
	
	void u_boundary(double *u_low, double *u_upp) override {
		u_low[0] = 0;
		u_upp[0] = 1;
	}
	
	void p_boundary(double *p_low, double *p_upp) override {
		
	}

	void var_boundary(double *x_low, double *x_upp) override {
		x_low[x_index(0,0)] = x_upp[x_index(0,0)] = 4;
		x_low[x_index(0,1)] = x_upp[x_index(0,1)] = 0;
	}

};


/////////////////////////////////////////////////////////////////////////////

int main(int argv, char* argc[]) {

	TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv,argc);
	TWfolder folder(&twparameter,0);

	FischenPhase ph(twparameter.NDIS);
	folder.Add(&ph);

	Viewer *viewer = nullptr;
	if (twparameter.PLOT) viewer = new Viewer(&twparameter);
	
	folder.Init();
	folder.Init(viewer);
	
	folder.Loop();

	delete viewer;
	
	return 0;
}
