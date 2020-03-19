/*----------------------------------------------------------------
 *
 *  Example: Raketenwagen mit Gitteranpassung
 *
 *----------------------------------------------------------------*/

#include "TWfolder.h"
#include "TWproblem.h"
#include "TransWORHP.h"

class RaketenPhase : public tw::TransWorhpProblem {
public:

	RaketenPhase(const tw::TWdimension &TWdata) : TransWorhpProblem(TWdata) {
		freetime = true;
	}
	
	void selectWindows(tw::Viewer *viewer) override {

		viewer->AddStateView(0,"x_1");
		viewer->AddStateView(1,"x_2");
		
		viewer->AddControlView(0,"u");
	}

	void p_init(double *p) override {
		p[0] = 3.0;
	}

	double obj() override {
		return p(0);
	}


	bool obj_structure(tw::DiffStructure &s) override {
		s(0,p_index(0) );
		return true;

	}

	void ode(double *dx, double t, const double *x, const double *u, const double *p) override {
		dx[0] = x[1] * p[0];
		dx[1] = u[0] * p[0];
	}


	// Optional
	bool ode_structure(tw::DiffStructure &s) override {
		s(0,x_indexode(1));
		s(0,p_indexode(0));

		s(1,u_indexode(0));
		s(1,p_indexode(0));
		
		return true;
	}


	bool ode_diff(tw::DiffStructure& s, double t, const double* x, const double* u, const double* p) override {
		s(0,x_indexode(1)) = p[0];
		s(1,u_indexode(0)) = p[0];
		return true;
	}
	
	bool ode_diff_p(tw::DiffStructure& s, double t, const double *x, const double *u, const double *p, int index) override {
		s(0,p_indexode(0)) = x[1];
		s(1,p_indexode(0)) = u[0];
		return true;
	}
		
	void u_boundary(double *u_low, double *u_upp) override {
		u_low[0] = -1.0;
		u_upp[0] = +1.0;

	}

	void p_boundary(double *p_low, double *p_upp) override {
		p_low[0] = 1.0;
		p_upp[0] = 10.0;
	}

	void var_boundary(double *x_low, double *x_upp) override {
		x_low[x_index(0,0)] = x_upp[x_index(0,0)] = 4.0;
		x_low[x_index(0,1)] = x_upp[x_index(0,1)] = -1.0;
		x_low[x_index(n_dis-1,0)] = x_upp[x_index(n_dis-1,0)] = 0.0;
		x_low[x_index(n_dis-1,1)] = x_upp[x_index(n_dis-1,1)] = 0.0;
	}
};



/////////////////////////////////////////////////////////////////////////////


int main(int argv, char* argc[]) {

	tw::TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv, argc);

	tw::Viewer *viewer = nullptr;
	if (twparameter.PLOT) viewer = new tw::Viewer(&twparameter);

	tw::TWdimension TWdim;
	TWdim.ID = "Raketenwagen mesh-refinement";
	TWdim.n_dis = twparameter.NDIS;
	TWdim.n_ode = 2;
	TWdim.n_ctrl = 1;
	TWdim.n_param = 1;

	RaketenPhase ph(TWdim);
	ph.setSolver(&twparameter);

	tw::TWfolder folder(&twparameter, 0);
	folder.Add(&ph);
	folder.Init();
	folder.Init(viewer);
	
	folder.meshRef();

	delete viewer;

	return 0;
}
