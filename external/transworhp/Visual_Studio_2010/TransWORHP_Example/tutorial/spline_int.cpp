/*----------------------------------------------------------------
 *
 *  Tutorial: Splineproblem in Bolza-Form (Integral-Form)
 *
 *----------------------------------------------------------------*/

#ifdef WIN32
#include "windows.h"
#endif

#include "TransWORHP.h"

using namespace std;

class SplineInt : public TransWorhp {
public:

	SplineInt(int dis) : TransWorhp("Spline / Bolza-Form",dis,2,1,0,0,0,1) {
		
		lagrange_weight[0] = 1;
	}

	double obj() {

		return 0;
	}

	bool obj_structure(DiffStructure &s) {

		// bei false wäre es dicht. Hier zunächst keine direkten Abhängigkeiten.
		// Werden über integral_structure bereitgestellt.
		return true;
	}

	bool obj_diff(DiffStructure &s) {

		return true;
	}
		
	void integral(double *f, double t, const double *x, const double *u,
			 const double *p) {
			
		f[0] = u[0]*u[0];
	}

	bool integral_structure(DiffStructure &s) {

		s(0, u_indexode(0));
		return true;
	}
	
	bool integral_diff(DiffStructure &s, double t, const double *x, const double *u, const double *p) {
		
		s(0, u_indexode(0)) = 2*u[0];
		return true;
	}
	
	
	void ode(double *dx, double t, const double *x, const double *u, const double *p) {

		dx[0] = x[1];
		dx[1] = u[0];
	}

	bool ode_structure(DiffStructure &s) {

		s(0, x_index(0,1)); // dx[0] / dx[1]
		s(1, u_index(0,0)); // dx[1] / du[0]
		return true;
	}

	bool ode_diff(DiffStructure &s, double t, const double *x,
			      const double *u, const double *p) {

		s(0, x_index(0,1))= 1;
		s(1, u_index(0,0))= 1;
		return true;
	}

	void u_boundary(double *u_low, double *u_upp) {

		u_low[0] = -6;
		u_upp[0] = +6;
	}

	void x_boundary(double *x_low, double *x_upp) {

	}

	void var_boundary(double *x_low, double *x_upp) {

		x_low[x_index(0,0)] = x_upp[x_index(0,0)] = 0;
		x_low[x_index(0,1)] = x_upp[x_index(0,1)] = 1;
	
		x_low[x_index(n_dis-1,0)] = x_upp[x_index(n_dis-1,0)] = 0;
		x_low[x_index(n_dis-1,1)] = x_upp[x_index(n_dis-1,1)] = 1;
	}

	void u_init(double *u, int i, int dis) {

		u[0] = -6 + (12.*i)/dis;
	}

};

///////////////////////////////////////////////////////////////////

int main(int argv, char* argc[]) {

	TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv,argc);
	TWfolder folder(&twparameter,0);

	SplineInt ph(twparameter.NDIS);
	folder.Add(&ph);

	Viewer *viewer = 0;
	if (twparameter.PLOT) viewer = new Viewer(&twparameter);

	folder.Init();
	folder.Init(viewer);
	folder.Loop();

	delete viewer;

	return 0;
}

