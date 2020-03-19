/*----------------------------------------------------------------
 *
 *  Example: Oszillator
 *
 *----------------------------------------------------------------*/

#include "TWfolder.h"
#include "TWproblem.h"
#include "TransWORHP.h"

class OszillatorPhase : public tw::TransWorhpProblem {
public:

	OszillatorPhase(const tw::TWdimension &TWdim) : TransWorhpProblem(TWdim) {
		freetime = true;
	}


	void selectWindows(tw::Viewer *viewer) override {

		viewer->AddStateView(0,"x_1");
		viewer->AddStateView(1,"x_2");
		
		viewer->AddControlView(0,"u");
	}


	void p_init(double *p) override {
		p[0] = 5.2;
	}


	double obj() override {
		return p(0);
	}

	bool obj_structure(tw::DiffStructure &s) override {
		s(0,p_index(0));
		return true;
	}

	bool obj_diff(tw::DiffStructure &s) override {
		s(0,p_index(0)) = 1;
		return true;
	}


	void ode(double *dx, double t, const double *x, const double *u, const double *p) override {
		dx[0] = x[1] * p[0];
		dx[1] = (-x[0] + u[0]) * p[0];
	}

	bool ode_structure(tw::DiffStructure &s) override {
		s(0,x_indexode(1));
		s(0,p_indexode(0));

		s(1,x_indexode(0));
		s(1,u_indexode(0));
		s(1,p_indexode(0));
		return true;
	}

	bool ode_diff(tw::DiffStructure &s, double t, const double *x, const double *u, const double *p) override {
		s(0,x_indexode(1)) = p[0];
	
		s(1,x_indexode(0)) = -p[0];
		s(1,u_indexode(0)) = p[0];
		return true;
	}
	
	bool ode_diff_p(tw::DiffStructure& s, double t, const double *x, const double *u, const double *p, int index) override {
		s(0,p_indexode(0)) = x[1];
		s(1,p_indexode(0)) = -x[0] + u[0];
		return true;
	}


	void u_boundary(double *u_low, double *u_upp) override {
		u_low[0] = -1;
		u_upp[0] = +1;
	}

	void p_boundary(double *p_low, double *p_upp) override {
		p_low[0] = 1;
		p_upp[0] = 10;
	}

	void var_boundary(double *x_low, double *x_upp) override {
		x_low[x_index(0,0)] = x_upp[x_index(0,0)] = 3;
		x_low[x_index(0,1)] = x_upp[x_index(0,1)] = 1;
		x_low[x_index(n_dis-1,0)] = x_upp[x_index(n_dis-1,0)] = 0;
		x_low[x_index(n_dis-1,1)] = x_upp[x_index(n_dis-1,1)] = 0;
	}
};


/////////////////////////////////////////////////////////////////////////////

int main(int argv, char* argc[]) {

	tw::TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv,argc);
	tw::TWfolder folder(&twparameter,0);

	tw::TWdimension TWdim;
	TWdim.ID = "Harmonischer Oszillator";
	TWdim.n_dis = twparameter.NDIS;
	TWdim.n_ode = 2;
	TWdim.n_ctrl = 1;
	TWdim.n_param = 1;
	TWdim.setMultinodes(3);

	OszillatorPhase ph(TWdim);
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
