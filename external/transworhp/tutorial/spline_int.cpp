/*----------------------------------------------------------------
 *
 *  Tutorial: Splineproblem in Bolza-Form (Integral-Form)
 *
 *----------------------------------------------------------------*/

#include "TWfolder.h"
#include "TWproblem.h"
#include "TransWORHP.h"

class SplineInt : public tw::TransWorhpProblem {
public:

	SplineInt(const tw::TWdimension &TWdata) : TransWorhpProblem(TWdata) {}

	void localinit() override {
		solver->lagrange_weight[0] = 1;
	}

	double obj() override {

		return 0;
	}

	bool obj_structure(tw::DiffStructure &s) override {

		// bei false waere es dicht. Hier zunaechst keine direkten Abhaengigkeiten.
		// Werden ueber integral_structure bereitgestellt.
		return true;
	}

	bool obj_diff(tw::DiffStructure &s) override {

		return true;
	}
		
	void integral(double *f, double t, const double *x, const double *u,
			 const double *p) override {
			
		f[0] = u[0]*u[0];
	}

	bool integral_structure(tw::DiffStructure &s) override {

		s(0, u_indexode(0));
		return true;
	}
	
	bool integral_diff(tw::DiffStructure &s, double t, const double *x, const double *u, const double *p) override {
		
		s(0, u_indexode(0)) = 2*u[0];
		return true;
	}
	
	
	void ode(double *dx, double t, const double *x, const double *u, const double *p) override {

		dx[0] = x[1];
		dx[1] = u[0];
	}

	bool ode_structure(tw::DiffStructure &s) override {

		s(0, x_index(0,1)); // dx[0] / dx[1]
		s(1, u_index(0,0)); // dx[1] / du[0]
		return true;
	}

	bool ode_diff(tw::DiffStructure &s, double t, const double *x,
			      const double *u, const double *p) override {

		s(0, x_index(0,1))= 1;
		s(1, u_index(0,0))= 1;
		return true;
	}

	void u_boundary(double *u_low, double *u_upp) override {

		u_low[0] = -6;
		u_upp[0] = +6;
	}

	void x_boundary(double *x_low, double *x_upp) override {

	}

	void var_boundary(double *x_low, double *x_upp) override {

		x_low[x_index(0,0)] = x_upp[x_index(0,0)] = 0;
		x_low[x_index(0,1)] = x_upp[x_index(0,1)] = 1;
	
		x_low[x_index(n_dis-1,0)] = x_upp[x_index(n_dis-1,0)] = 0;
		x_low[x_index(n_dis-1,1)] = x_upp[x_index(n_dis-1,1)] = 1;
	}

	void u_init(double *u, int i, int dis) override {

		u[0] = -6 + (12.*i)/dis;
	}

};

///////////////////////////////////////////////////////////////////

int main(int argv, char* argc[]) {

	tw::TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv, argc);
	tw::TWfolder folder(&twparameter, 0);

	tw::TWdimension TWdim;
	TWdim.ID = "Spline / Bolza-Form";
	TWdim.n_ode = 2;
	TWdim.n_ctrl = 1;
	TWdim.n_integral = 1;
	TWdim.n_dis = twparameter.NDIS;

	SplineInt ph(TWdim);
	folder.Add(&ph);

	tw::Viewer *viewer = nullptr;
	if (twparameter.PLOT) viewer = new tw::Viewer(&twparameter);

	folder.Init();
	folder.Init(viewer);
	folder.Loop();

	delete viewer;

	return 0;
}
