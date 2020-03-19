/*----------------------------------------------------------------
 *
 *  Example: Eismann
 *
 *----------------------------------------------------------------*/

#ifdef WIN32
#include "windows.h"
#endif

#include "TransWORHP.h"

using namespace std;

double preis(double t);

class EismannPhase : public TransWorhp {
public:

	EismannPhase(int dis) : TransWorhp("Eismann", dis,2,1,1,0,0) {
		freetime = 1;
	}

	string GetXTitle(int d) override {
		if (d==0) return "Kapital";
		if (d==1) return "Eisbestand";
	}

	string GetUTitle(int d) override {
		if (d==0) return "Steuerung";
	}

	double preis(double t) {
		if (t >= 0 && t < 4) {
			return 6.0 + 0.5*t;
		} else if (t >= 4 && t < 6) {
			return 4.0 + t;
		} else {
			return 10.0;
		}
	}
	
	
	void p_init(double *p) override {
		p[0] = 8;
	}

	void u_init(double *u, int i, int dis) override {
		if (i < ((16.0/3)/8)*dis)
			u[0] = 1.0;
		else
			u[0] = -1.0;
	}

	void x_init(double *x, int i, int dis) override {
		x[0] = 10.0;
		x[1] = 2.0;
	}
	
	double obj() override {
		return -x(n_dis-1,0) - preis(p(0))*x(n_dis-1,1);
	}

	bool obj_structure(DiffStructure &s) override {
		//s(0,p_index(0) );
		return false;
	}
	
	bool obj_diff(DiffStructure &s) override {
		//s(0,p_index(0) ) = 1;
		return false;
	}
	

	void ode(double *dx, double t, const double *x, const double *u, const double *p) override {
		dx[0] = (-0.25*x[1]-preis(t*p[0])*u[0]) * p[0];
		dx[1] = u[0] * p[0];
	}

	bool ode_structure(DiffStructure &s) override {
		//s(0,x_indexode(1));
		//s(0,p_indexode(0));

		//s(1,u_indexode(0));
		//s(1,p_indexode(0));
		return false;
	}

	bool ode_diff(DiffStructure& s, double t, const double* x, const double* u, const double* p) override {
		//s(0,x_indexode(1)) = p[0];
		//s(1,u_indexode(0)) = p[0];
		return false;
	}
	
	bool ode_diff_p(DiffStructure& s, double t, const double *x, const double *u, const double *p, int index) override {
		//s(0,p_indexode(0)) = x[1];
		//s(1,p_indexode(0)) = u[0];
		return false;
	}
		
		
	void u_boundary(double *u_low, double *u_upp) override {
		u_low[0] = -1;
		u_upp[0] = +1;
	}

	void p_boundary(double *p_low, double *p_upp) override {
		p_low[0] = 8;
		p_upp[0] = 8;
	}

	void var_boundary(double *x_low, double *x_upp) override {
		x_low[x_index(0,0)] = x_upp[x_index(0,0)] = 10;
		x_low[x_index(0,1)] = x_upp[x_index(0,1)] = 2;
	}

};


/////////////////////////////////////////////////////////////////////////////

int main(int argv, char* argc[]) {

	TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv,argc);
	TWfolder folder(&twparameter,0);

	EismannPhase ph(twparameter.NDIS);
	folder.Add(&ph);

	Viewer *viewer = nullptr;
	if (twparameter.PLOT) viewer = new Viewer(&twparameter);
	
	folder.Init();
	folder.Init(viewer);
	
	folder.Loop();

	delete viewer;
	
	return 0;
}
