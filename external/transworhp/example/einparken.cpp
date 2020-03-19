/*----------------------------------------------------------------
 *
 *  Example: Einparken
 *
 *----------------------------------------------------------------*/

#include "TWfolder.h"
#include "TWproblem.h"
#include "TransWORHP.h"

double phasenPlot (int &len, int &ndgl, int &nsteuer, double *t,double *x, double *u, int &i, int &index);

class EinparkenPhase : public tw::TransWorhpProblem {
public:

	EinparkenPhase(int dis) : TransWorhpProblem(tw::TWdimension("Einparken", dis,6,2,1,0,0)) {
		freetime = true;
	}

	void OpenWindows(tw::Viewer *viewer) override {
		viewer->PhasePlot("Phasen-Plot", phasenPlot,0,1);
	}

	void selectWindows(tw::Viewer *viewer) override {

		viewer->AddStateView(0,"x");
		viewer->AddStateView(1,"y");
		viewer->AddStateView(2,"v");
		viewer->AddStateView(3,"theta");
		viewer->AddStateView(4,"phi");
		viewer->AddStateView(5,"hilfe");
		
		viewer->AddControlView(0,"Beschleunigung");
		viewer->AddControlView(1,"sigma");
	}
	
	void p_init(double *p) override {
		p[0] = 2.0;
	}
	
	double obj() override {
		return p(0) + 0.0*x(n_dis-1, 5);
	}

	bool obj_structure(tw::DiffStructure &s) override {
		s(0,x_index(n_dis-1,5));
		s(0,p_index(0));
		return true;
	}
	
	bool obj_diff(tw::DiffStructure &s) override {
		s(0,x_index(n_dis-1,5)) = 0.0;
		s(0,p_index(0)) = 1;
		return true;
	}
	
	void ode(double *dx, double t, const double *x, const double *u, const double *p) override {
		dx[0] = x[2] * cos(x[3]) * p[0];
		dx[1] = x[2] * sin(x[3]) * p[0];
		dx[2] = u[0] * p[0];
		dx[3] = x[2] * sin(x[4]) * p[0] / 1; // b = 1
		dx[4] = u[1] * p[0];
		dx[5] = u[0]*u[0] + u[1]*u[1];
	}

	bool ode_structure(tw::DiffStructure &s) override {
		s(0,x_indexode(2));
		s(0,x_indexode(3));
		s(0,p_indexode(0));

		s(1,x_indexode(2));
		s(1,x_indexode(3));
		s(1,p_indexode(0));

		s(2,u_indexode(0));
		s(2,p_indexode(0));

		s(3,x_indexode(2));
		s(3,x_indexode(4));
		s(3,p_indexode(0));

		s(4,u_indexode(1));
		s(4,p_indexode(0));

		s(5,u_indexode(0));
		s(5,u_indexode(1));
		return true;
	}

	bool ode_diff(tw::DiffStructure& s, double t, const double* x, const double* u, const double* p) override {
		s(0,x_indexode(2)) = cos(x[3]) * p[0];
		s(0,x_indexode(3)) = -x[2] * sin(x[3]) * p[0];

		s(1,x_indexode(2)) = sin(x[3]) * p[0];
		s(1,x_indexode(3)) = x[2] * cos(x[3]) * p[0];		

		s(2,u_indexode(0)) = p[0];

		s(3,x_indexode(2)) = sin(x[4]) * p[0] / 1; // b = 1
		s(3,x_indexode(4)) = x[2] * cos(x[4]) * p[0] / 1; // b = 1

		s(4,u_indexode(1)) = p[0];

		s(5,u_indexode(0)) = 2*u[0];
		s(5,u_indexode(1)) = 2*u[1];
		return true;
	}
	
	bool ode_diff_p(tw::DiffStructure& s, double t, const double *x, const double *u, const double *p, int index) override {
		s(0,p_indexode(0)) = x[2] * cos(x[3]);
		s(1,p_indexode(0)) = x[2] * sin(x[3]);
		s(2,p_indexode(0)) = u[0];
		s(3,p_indexode(0)) = x[2] * sin(x[4]) / 1; // b = 1
		s(4,p_indexode(0)) = u[1];
		return true;
	}
	
	void x_boundary(double *x_low, double *x_upp) override {
		x_low[2] = 0.0;
		x_upp[2] = 10.0;
		x_low[4] = -M_PI/2; // phi_max
		x_upp[4] = M_PI/2;
		x_low[5] = 0.0;
	}
	
	void u_boundary(double *u_low, double *u_upp) override {
		u_low[0] = -10.0;
		u_upp[0] = 10.0;
		u_low[1] = -M_PI/2; // sigma_max
		u_upp[1] = M_PI/2;
	}

	void p_boundary(double *p_low, double *p_upp) override {
		p_low[0] = 1.0;
		p_upp[0] = 10.0;
	}

	void var_boundary(double *x_low, double *x_upp) override {
		x_low[x_index(0,0)] = x_upp[x_index(0,0)] = 0.0;
		x_low[x_index(0,1)] = x_upp[x_index(0,1)] = 0.0;
		x_low[x_index(0,2)] = x_upp[x_index(0,2)] = 0.0;
		x_low[x_index(0,3)] = x_upp[x_index(0,3)] = 0.0;
		x_low[x_index(0,4)] = x_upp[x_index(0,4)] = 0.0;
		x_low[x_index(n_dis-1,0)] = x_upp[x_index(n_dis-1,0)] = 5.0;
		x_low[x_index(n_dis-1,1)] = x_upp[x_index(n_dis-1,1)] = -2.0;
		x_low[x_index(n_dis-1,2)] = x_upp[x_index(n_dis-1,2)] = 0.0;
		x_low[x_index(n_dis-1,3)] = x_upp[x_index(n_dis-1,3)] = -M_PI/2;
		x_low[x_index(n_dis-1,4)] = x_upp[x_index(n_dis-1,4)] = 0.0;
	}

};

double phasenPlot (int &len, int &ndgl, int &nsteuer, double *t,double *x, double *u, int &i, int &index) {
	if (index==0) {			return x[i*(ndgl+nsteuer)+1];
	} else if (index==1) {	return x[i*(ndgl+nsteuer)];
	} else {				return 0.0;}
}

/////////////////////////////////////////////////////////////////////////////

int main(int argv, char* argc[]) {

	tw::TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv,argc);
	tw::TWfolder folder(&twparameter,0);

	EinparkenPhase ph(twparameter.NDIS);
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
