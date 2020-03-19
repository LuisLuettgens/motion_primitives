/*----------------------------------------------------------------
 *
 * Example: Chemischer Reaktor
 *
 *----------------------------------------------------------------*/

#include "TWfolder.h"
#include "TWproblem.h"
#include "TransWORHP.h"

class ChemiePhase : public tw::TransWorhpProblem {
public:

	ChemiePhase(const tw::TWdimension &TWdata) : TransWorhpProblem(TWdata) {}

	void selectWindows(tw::Viewer *viewer) override {

		viewer->AddStateView(0,"x_1");
		viewer->AddStateView(1,"x_2");
		
		viewer->AddControlView(0,"u");
	}
	
	double obj() override {
		return -x(n_dis-1,1);
	}

	bool obj_structure(tw::DiffStructure &s) override {
		s(0,x_index(n_dis-1,1));
		return true;
	}

	bool obj_diff(tw::DiffStructure &s) override {
		s(0,x_index(n_dis-1,1)) = -1.0;
		return true;
	}
	
	
	void ode(double *dx, double t, const double *x, const double *u, const double *p) override {
		dx[0] = -u[0]*x[0] +   u[0]*u[0]*x[1];
		dx[1] =  u[0]*x[0] - 3*u[0]*u[0]*x[1];
	}

	bool ode_structure(tw::DiffStructure &s) override {
		s(0,u_indexode(0));
		s(0,x_indexode(0));
		s(0,x_indexode(1));

		s(1,u_indexode(0));
		s(1,x_indexode(0));
		s(1,x_indexode(1));
		return true;
	}

	bool ode_diff(tw::DiffStructure& s, double t, const double* x, const double* u, const double* p) override {
		s(0,u_indexode(0)) = -x[0] + 2*u[0]*x[1];
		s(0,x_indexode(0)) = -u[0];
		s(0,x_indexode(1)) = u[0]*u[0];

		s(1,u_indexode(0)) = x[0] - 6*u[0]*x[1];
		s(1,x_indexode(0)) = u[0];
		s(1,x_indexode(1)) = -3*u[0]*u[0];
		return true;
	}


	void u_boundary(double *u_low, double *u_upp) override {
		u_low[0] =  0.0;
		u_upp[0] = +1.0;
	}

	void var_boundary(double *x_low, double *x_upp) override {
		x_upp[x_index(0,0)] = x_low[x_index(0,0)] = 1.0;
		x_upp[x_index(0,1)] = x_low[x_index(0,1)] = 0.0;
	}

};

/////////////////////////////////////////////////////////////////////////////

int main(int argv, char* argc[]) {

	tw::TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv, argc);
	tw::TWfolder folder(&twparameter, 0);

	tw::TWdimension TWdim;
	TWdim.ID = "Chemiereaktor";
	TWdim.n_dis = twparameter.NDIS;
	TWdim.n_ode = 2;
	TWdim.n_ctrl = 1;

	ChemiePhase ph(TWdim);
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
