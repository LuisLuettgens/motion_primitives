/*----------------------------------------------------------------
 *
 *  Example: Eismann
 *
 *----------------------------------------------------------------*/

#include "TWfolder.h"
#include "TWproblem.h"
#include "TransWORHP.h"

class EismannPhase : public tw::TransWorhpProblem {
public:

	EismannPhase(int dis) : TransWorhpProblem(tw::TWdimension("Eismann", dis,2,1,1,0,0)) {
		freetime = true;
	}
	
	void selectWindows(tw::Viewer *viewer) override {

		viewer->AddStateView(0,"Kapital");
		viewer->AddStateView(1,"Eisbestand");
		
		viewer->AddControlView(0,"Steuerung");
	}
	
	double preis(double t) const {
		if (t >= 0 && t < 4) {
			return 6.0 + 0.5*t;
		} else if (t >= 4 && t < 6) {
			return 4.0 + t;
		} else {
			return 10.0;
		}
	}
	
	void p_init(double *p) override {
		p[0] = 8.0;
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

	bool obj_structure(tw::DiffStructure &s) override {
		s(0,x_index(n_dis-1,0));
		s(0,x_index(n_dis-1,1));
		s(0,p_index(0));
		return true;
	}
	
	bool obj_diff(tw::DiffStructure &s) override {
		//s(0,p_index(0) ) = 1;
		return false;
	}
	

	void ode(double *dx, double t, const double *x, const double *u, const double *p) override {
		dx[0] = (-0.25*x[1]-preis(t*p[0])*u[0]) * p[0];
		dx[1] = u[0] * p[0];
	}

	bool ode_structure(tw::DiffStructure &s) override {
		//s(0,x_indexode(1));
		//s(0,p_indexode(0));

		//s(1,u_indexode(0));
		//s(1,p_indexode(0));
		return false;
	}

	bool ode_diff(tw::DiffStructure& s, double t, const double* x, const double* u, const double* p) override {
		//s(0,x_indexode(1)) = p[0];
		//s(1,u_indexode(0)) = p[0];
		return false;
	}
	
	bool ode_diff_p(tw::DiffStructure& s, double t, const double *x, const double *u, const double *p, int index) override {
		//s(0,p_indexode(0)) = x[1];
		//s(1,p_indexode(0)) = u[0];
		return false;
	}
		
		
	void u_boundary(double *u_low, double *u_upp) override {
		u_low[0] = -1.0;
		u_upp[0] = +1.0;
	}

	void p_boundary(double *p_low, double *p_upp) override {
		p_low[0] = 8.0;
		p_upp[0] = 8.0;
	}

	void var_boundary(double *x_low, double *x_upp) override {
		x_low[x_index(0,0)] = x_upp[x_index(0,0)] = 10.0;
		x_low[x_index(0,1)] = x_upp[x_index(0,1)] = 2.0;
	}

};

/////////////////////////////////////////////////////////////////////////////

int main(int argv, char* argc[]) {

	tw::TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv,argc);
	tw::TWfolder folder(&twparameter,0);

	EismannPhase ph(twparameter.NDIS);
	ph.setSolver(&twparameter);
	folder.Add(&ph);

	tw::Viewer *viewer = nullptr;
	if (twparameter.PLOT) viewer = new tw::Viewer(&twparameter);
	
	folder.Init();
	folder.Init(viewer);
	
	folder.Loop();

	delete viewer;
	
	return 0;
}
