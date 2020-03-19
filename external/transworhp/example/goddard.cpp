/*----------------------------------------------------------------
 *
 *  Example: Goddard Problem / Hoehenrakete
 *
 *----------------------------------------------------------------*/

#include "TWfolder.h"
#include "TWproblem.h"
#include "TransWORHP.h"

class Goddard : public tw::TransWorhpProblem {
public:

	Goddard(int dis) : TransWorhpProblem(tw::TWdimension("Hoehenrakete", dis, 3, 1, 1, 0, 0)) {
		freetime = true;
	}

	void selectWindows(tw::Viewer *viewer) override {

		viewer->AddStateView(0,"h - Hoehe");
		viewer->AddStateView(1,"v - Geschwindigkeit");
		viewer->AddStateView(2,"m - Masse");
		
		viewer->AddControlView(0,"u - Steuerung");
	}
	
        void x_init(double *x, int i, int dis) override {
		x[0] = 1.0;
		x[1] = 1.0;
		x[2] = 1.0;
	}

	void p_init(double *p) override {
		p[0] = 211.6;
	}

	void u_init(double *u, int i, int dis) override {
		u[0] = 1.0;
	}

	double obj() override {
		return -x(n_dis-1, 0);
	}

	bool obj_structure(tw::DiffStructure &s) override {
		s(0, x_index(n_dis-1, 0));
		return true;
	}

	bool obj_diff(tw::DiffStructure &s) override {
		s(0, x_index(n_dis-1, 0) ) = -1.0;
		return true;
	}

	void ode(double *dx, double t, const double *x, const double *u, const double *p) override {
		
		const double c = 2060.103;
		const double alpha = 1.22710e-2;
		const double beta = 1.45010e-4;
		const double g0 = 9.810;
		const double r0 = 6.37110e6;
		
		double g = g0 * (r0/(r0+x[0]))*(r0/(r0+x[0]));
		double d = alpha*x[1]*x[1] * exp(-beta*x[0]);
	  
		dx[0] = x[1];
		dx[1] = 1.0/x[2]*(u[0]*c - d) - g;
		dx[2] = -u[0];

		// free final time
		dx[0] *= p[0];
		dx[1] *= p[0];
		dx[2] *= p[0];
	}

	bool ode_structure(tw::DiffStructure &s) override {
		/*
		s(0, x_indexode(1));

		s(1, x_indexode(0));
		s(1, x_indexode(1));
		s(1, x_indexode(2));
		s(1, u_indexode(0));

		s(2, u_indexode(0));
		
		s(0, p_indexode(0));
		s(1, p_indexode(0));
		s(2, p_indexode(0));
		*/
		return false;
	}

	bool ode_diff(tw::DiffStructure& s, double t, const double* x, const double* u, const double* p) override {
		
		/*
		s(0, x_indexode(1)) = 1.0;

		s(1, x_indexode(0));
		s(1, x_indexode(1));
		s(1, x_indexode(2));
		s(1, u_indexode(0));

		s(2, u_indexode(0)) = -1.0;
		
		s(0, p_indexode(0)) = x[1];
		s(1, p_indexode(0)) = 1.0/x[2]*(u[0]*c - d) - g;
		s(2, p_indexode(0)) = -u[0];
		*/
		
		return false;
	}

	bool ode_diff_p(tw::DiffStructure& s, double t, const double *x, const double *u, const double *p, int index) override {
		
		return false;
	}

	void x_boundary(double *x_low, double *x_upp) override {
	  
	}

	void u_boundary(double *u_low, double *u_upp) override {
		u_low[0] = 0.0;
		u_upp[0] = 10.0;
	}

	void p_boundary(double *p_low, double *p_upp) override {
		p_low[0] = 200.0;
		p_upp[0] = 300.0;
	}

	void var_boundary(double *x_low, double *x_upp) override {

		x_low[x_index(0,0)] = x_upp[x_index(0,0)] = 0.0;
		x_low[x_index(0,1)] = x_upp[x_index(0,1)] = 0.0;
		x_low[x_index(0,2)] = x_upp[x_index(0,2)] = 214.839;

		x_low[x_index(n_dis-1,2)] = x_upp[x_index(n_dis-1,2)] = 67.9833;
	}
	
};

int main(int argv, char* argc[]) {

	tw::TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv, argc);

	tw::Viewer *viewer = nullptr;
	if (twparameter.PLOT) viewer = new tw::Viewer(&twparameter);

	Goddard ph(twparameter.NDIS);
	ph.setSolver(&twparameter);

	tw::TWfolder folder(&twparameter, 0);
	folder.Add(&ph);
	folder.Init();
	folder.Init(viewer);
	
	folder.meshRef();
	
	delete viewer;

	return 0;
}
