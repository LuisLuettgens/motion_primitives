/*----------------------------------------------------------------
 *
 *  Tutorial: Splineproblem
 *
 *----------------------------------------------------------------*/

#ifdef WIN32
#include "windows.h"
#endif

#include "TransWORHP.h"

using namespace std;

class SplineTest : public TransWorhp {
public:

	SplineTest(int dis) : TransWorhp("Spline",dis,3,1,0,0,1) {}

	double obj() {

		return x(n_dis-1,2);
	}

	bool obj_structure(DiffStructure &s) {

		s(0, x_index(n_dis-1,2));
		return true;
	}

	bool obj_diff(DiffStructure &s) {

		s(0, x_index(n_dis-1,2)) = 1;
		return true;
	}

	void ode(double *dx, double t, const double *x, const double *u,
			 const double *p) {

		dx[0] = x[1];
		dx[1] = u[0];
		dx[2] = u[0]*u[0];
	}

	bool ode_structure(DiffStructure &s) {

		s(0, x_indexode(1)); // dx[0] / dx[1]
		s(1, u_indexode(0)); // dx[1] / du[0]
		s(2, u_indexode(0)); // dx[2] / du[0]
		return true;
	}

	bool ode_diff(DiffStructure &s, double t, const double *x,
			      const double *u, const double *p) {

		s(0, x_indexode(1))= 1;
		s(1, u_indexode(0))= 1;
		s(2, u_indexode(0))= 2*u[0];
		return true;
	}

	void u_boundary(double *u_low, double *u_upp) {

		u_low[0] = -6;
		u_upp[0] = +6;
	}

	void x_boundary(double *x_low, double *x_upp) {

		x_low[2] = 0;
	}


	void var_boundary(double *x_low, double *x_upp) {

		x_low[x_index(0,0)] = x_upp[x_index(0,0)] = 0;
		x_low[x_index(0,1)] = x_upp[x_index(0,1)] = 1;
		x_low[x_index(0,2)] = x_upp[x_index(0,2)] = 0;

		x_low[x_index(n_dis-1,0)] = x_upp[x_index(n_dis-1,0)] = 0;
		x_low[x_index(n_dis-1,1)] = x_upp[x_index(n_dis-1,1)] = 1;
	}

/*
	void rand(double *r) {

		r[0] = x(0,0);
		r[1] = x(0,1) - 1;
		r[2] = x(0,2);

		r[3] = x(n_dis-1,0);
		r[4] = x(n_dis-1,1) - 1;
		
		r[5] = x__(3,0)-x__(5,0);
		r[6] = x__(4,0)-x__(8,0);
	}

	bool rand_structure(DiffStructure &s) {

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
	
	bool rand_diff(DiffStructure &s) {
	
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
	}*/

	
	void neben(double *c, double t, const double *x, 
			   const double *u, const double *p) {
		
		c[0] = x[0] + x[1];	
		
		//t???
	}
	
    void neben_boundary(double *c_low, double *c_upp) {
		
		c_low[0] = -.3;
		c_upp[0] = 1;
		
	}

	bool neben_structure(DiffStructure &s) {

		s(0,x_index(0,0));
    	s(0,x_index(0,1));
        return true;
    }
    
    bool neben_diff(DiffStructure &s, double t, const double *x,
					const double *u, const double *p) {
	
		s(0,x_index(0,0)) = 1;
		s(0,x_index(0,1)) = 1;
		return true;
	}
	

	void u_init(double *u, int i, int dis) {

		u[0] = -6 + (12.*i)/dis;
	}

	void terminate() {
		
/*		double SplineKoeff[][];
		
		for (int k=0; k<n_dis-1; k++) { // für jedes Intervall
		
			for (int i=0; i<n_ode; i++) { // Spline x_i anlegen
				
				double x_k[i] = x(k,i);
				double x_k1 = x(k+1,i);
				double u_k = u(k,i);
				double u_k1 = u(k+1,i);
				
				double f_k = new double[n_ode];
				//ode(double *dx, double t, const double *x, const double *u, const double *p)
				ode(f_k, T[k], x_k, u_k, p);

			}
			
			for (int i=0; i<n_ctrl; i++) { // Spline u_i anlegen
				
			}
		}
*/		
		
	}
};

///////////////////////////////////////////////////////////////////

int main(int argv, char* argc[]) {

	TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv,argc);
	TWfolder folder(&twparameter,0);

	SplineTest ph(twparameter.NDIS);
	folder.Add(&ph);

	Viewer *viewer = 0;
	if (twparameter.PLOT) viewer = new Viewer(&twparameter);

	folder.Init();
	folder.Init(viewer);
	folder.Loop();

	delete viewer;

	return 0;
}

