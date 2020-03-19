/*----------------------------------------------------------------
 *
 *  Example: Splineproblem mit Beschr�nkung
 *
 *----------------------------------------------------------------*/

#ifdef WIN32
#include "windows.h"
#endif

#include "TransWORHP.h"
#include "Viewer.h"
#include "TWfolder.h"

using namespace std;

class SplinePhase : public TransWorhp {
public:

	SplinePhase(int dis) : TransWorhp("Spline",dis,3,1,0,0,0) {}

	void init() {

		for (int i=0;i<n_dis;i++) {
			X[u_index(i,0)] = .5;
			//X[x_index(i,1)] =  + .3 * i/n_dis;
		}

	}

	double obj() {
		return x(n_dis-1,2);
	}

	bool obj_structure(DiffStructure &s) {
		s(0, x_index(n_dis-1,2));
		return true;
	}

	void ode(double *dx, double t, const double *x, const double *u, const double *p) {

// void ode(TransWORHP::Vektor &dx, double t, const TransWORHP::Vektor &x, const TransWORHP::Vektor &u,
//  const TransWORHP::Vektor &p) {

		dx[0] = x[1];
		dx[1] = u[0];
		dx[2] = u[0]*u[0];

		//cout << x ;
		//cout << u;
		//cout << dx;
	}

	// Optional
	bool ode_structure(DiffStructure &s) {

		s(0,x_index(0,1)); // dx[0] / dx[1]
		s(1,u_index(0,0)); // dx[1] / du[0]
		s(2,u_index(0,0)); // dx[2] / du[0]


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
		x_low[x_index(0,2)] = x_upp[x_index(0,2)] = 0;

		x_low[x_index(n_dis-1,0)] = x_upp[x_index(n_dis-1,0)] = 0;
		x_low[x_index(n_dis-1,1)] = x_upp[x_index(n_dis-1,1)] = 1;

	}


	/*void rand(double *r) {
		r[0] = x(0,0);//*x(0,0)*x(n_dis-1,0)*x(n_dis-1,0);
		r[1] = x(0,1) - 1;
		r[2] = x(0,2);

		r[3] = x(n_dis-1,0);
		r[4] = x(n_dis-1,1) - 1;

		// r_boundary
	}

	bool rand_structure(DiffStructure &s) {

		s(0,x_index(0,0)); //
		s(1,x_index(0,1)); //
		s(2,x_index(0,2)); //
		s(3,x_index(n_dis-1,0)); //
		s(4,x_index(n_dis-1,1)); //

		return true;
	}*/

};


class ConstrSplinePhase : public SplinePhase {
public:

	ConstrSplinePhase(int dis) : SplinePhase(dis) {
		id="Constrained Spline";
	}

	void x_boundary(double *x_low, double *x_upp) {

		SplinePhase::x_boundary(x_low,x_upp);
		x_upp[0] = 0.1;
	}

	void var_boundary(double *x_low, double *x_upp) {

		SplinePhase::var_boundary(x_low,x_upp);
		x_low[x_index(0,0)] = x_upp[x_index(0,0)] = 0;
		x_low[x_index(0,1)] = x_upp[x_index(0,1)] = 1;
		x_low[x_index(0,2)] = x_upp[x_index(0,2)] = 0;

		x_low[x_index(n_dis-1,0)] = x_upp[x_index(n_dis-1,0)] = 0;
		x_low[x_index(n_dis-1,1)] = x_upp[x_index(n_dis-1,1)] = -1;

	}
	/*void rand(double *r) {

		SplinePhase::rand(r);
		r[4] = x(n_dis-1,1) + 1;
	}*/


};

/////////////////////////////////////////////////////////////////////////////

int main(int argv, char* argc[]) {

	TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv,argc);
	TWfolder folder(&twparameter,0);

	ConstrSplinePhase ph(twparameter.NDIS);
	folder.Add(&ph);

	Viewer *viewer = 0;
	if (twparameter.PLOT) viewer = new Viewer(&twparameter);

	folder.Init();
	folder.Init(viewer);
	folder.Loop();

	delete viewer;

	return 0;
}

