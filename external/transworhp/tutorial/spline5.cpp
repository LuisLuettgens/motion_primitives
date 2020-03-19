/*----------------------------------------------------------------
 *
 *  Tutorial: Splineproblem mit Angabe der Struktur und Ableitung
 *            + Neben- und Randbedingungen
 *
 *----------------------------------------------------------------*/

#include "TWfolder.h"
#include "TWproblem.h"
#include "TransWORHP.h"

class Spline5 : public tw::TransWorhpProblem {
public:

	Spline5(const tw::TWdimension &TWdata) : TransWorhpProblem(TWdata) {}

	double obj() override {

		return x(n_dis-1,2);
	}

	bool obj_structure(tw::DiffStructure &s) override {

		s(0, x_index(n_dis-1,2));
		return true;
	}

	bool obj_diff(tw::DiffStructure &s) override {

		s(0, x_index(n_dis-1,2)) = 1;
		return true;
	}

	void ode(double *dx, double t, const double *x, const double *u,
			 const double *p) override {

		dx[0] = x[1];
		dx[1] = u[0];
		dx[2] = u[0]*u[0];
	}

	bool ode_structure(tw::DiffStructure &s) override {

		s(0, x_indexode(1)); // dx[0] / dx[1]
		s(1, u_indexode(0)); // dx[1] / du[0]
		s(2, u_indexode(0)); // dx[2] / du[0]
		return true;
	}

	bool ode_diff(tw::DiffStructure &s, double t, const double *x,
			      const double *u, const double *p) override {

		s(0, x_indexode(1))= 1;
		s(1, u_indexode(0))= 1;
		s(2, u_indexode(0))= 2*u[0];
		return true;
	}

	void u_boundary(double *u_low, double *u_upp) override {

		u_low[0] = -6;
		u_upp[0] = +6;
	}

	void x_boundary(double *x_low, double *x_upp) override {

		x_low[2] = 0;
	}

/*
	void var_boundary(double *x_low, double *x_upp) {

		x_low[x_index(0,0)] = x_upp[x_index(0,0)] = 0;
		x_low[x_index(0,1)] = x_upp[x_index(0,1)] = 1;
		x_low[x_index(0,2)] = x_upp[x_index(0,2)] = 0;

		x_low[x_index(n_dis-1,0)] = x_upp[x_index(n_dis-1,0)] = 0;
		x_low[x_index(n_dis-1,1)] = x_upp[x_index(n_dis-1,1)] = 1;
	}
*/

	void rand(double *r) override {

		r[0] = x(0,0);
		r[1] = x(0,1) - 1;
		r[2] = x(0,2);

		r[3] = x(n_dis-1,0);
		r[4] = x(n_dis-1,1) - 1;
		
		r[5] = x__(3,0)-x__(5,0);
		r[6] = x__(4,0)-x__(8,0);
	}

	bool rand_structure(tw::DiffStructure &s) override {

		s(0,x_index(0,0));
		s(1,x_index(0,1));
		s(2,x_index(0,2));
		s(3,x_index(n_dis-1,0));
		s(4,x_index(n_dis-1,1));
		
		s(5,x_index__(3,0));
		s(5,x_index__(5,0));
		
		s(6,x_index__(4,0));
		s(6,x_index__(8,0));
		
		return true;
	}
	
	bool rand_diff(tw::DiffStructure &s) override {
	
		s(0,x_index(0,0)) = 1;
		s(1,x_index(0,1)) = 1;
		s(2,x_index(0,2)) = 1;
		s(3,x_index(n_dis-1,0)) = 1;
		s(4,x_index(n_dis-1,1)) = 1;
	
		s(5,x_index__(3,0)) = 1;
		s(5,x_index__(5,0)) = -1;
		
		s(6,x_index__(4,0)) = 1;
		s(6,x_index__(8,0)) = -1;;
	
		return true;
	}

	
	void neben(double *c, double t, const double *x, 
			   const double *u, const double *p) override {
		
		c[0] = x[0] + x[1];	
	}
	
    	void neben_boundary(double *c_low, double *c_upp) override {
		
		c_low[0] = -.4;
		c_upp[0] = 1;
	}

	bool neben_structure(tw::DiffStructure &s) override {

		s(0,x_index(0,0));
    		s(0,x_index(0,1));
        	return true;
    	}
    
    bool neben_diff(tw::DiffStructure &s, double t, const double *x,
					const double *u, const double *p) override {
	
		s(0,x_index(0,0)) = 1;
		s(0,x_index(0,1)) = 1;
		return true;
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
	TWdim.ID = "Spline5";
	TWdim.n_ode = 3;
	TWdim.n_ctrl = 1;
	TWdim.n_rand = 7;
	TWdim.n_neben = 1;
	TWdim.n_dis = twparameter.NDIS;

	Spline5 ph(TWdim);
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
