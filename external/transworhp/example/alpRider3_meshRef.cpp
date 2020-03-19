/*----------------------------------------------------------------
*
* Alp Rider aus Betts
*
*----------------------------------------------------------------*/

#include "TWfolder.h"
#include "TWproblem.h"
#include "TransWORHP.h"

class AlpRiderPhase : public tw::TransWorhpProblem {
public:

	AlpRiderPhase(int dis) : TransWorhpProblem(tw::TWdimension("Alp Rider", dis, 4, 2, 0, 0, 1, 1)) {
		//freetime = true;
	}
	
	void localinit() {
		solver->lagrange_weight[0] = 1.0;
	}

	double peak(double t, double a, double b) {
		return exp(-b*(t - a)*(t - a));
	}

	void x_init(double *x, int /*i*/, int /*dis*/) {
		x[0] = 2;
		x[1] = 1;
		x[2] = 1;
		x[3] = 1;
	}

	void u_init(double *u, int /*i*/, int /*dis*/) {
		u[0] = 1;
		u[1] = 1;
	}

	double obj() {
		return 0;//x(n_dis-1,4);
	}

	bool obj_structure(tw::DiffStructure &/*s*/) {
		//	s(0,x_index(n_dis-1,4) );
		return true;

	}
	/*
	bool obj_diff(DiffStructure &s) {
	s(0,x_index(n_dis-1,8)) = 1 ;
	s(0,p_index(0)) = .1 ;
	return true;
	}
	*/

	void integral(double *f, double /*t*/, const double *x, const double *u, const double */*p*/) {

		const double c = 20.0;

		f[0] = ((100.0 * (x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3]) +
			0.01 * (u[0] * u[0] + u[1] * u[1]))) * c;

	}

	bool integral_structure(tw::DiffStructure &s) {

		s(0, x_indexode(0));
		s(0, x_indexode(1));
		s(0, x_indexode(2));
		s(0, x_indexode(3));
		s(0, u_indexode(0));
		s(0, u_indexode(1));
		return true;
	}

	bool integral_diff(tw::DiffStructure &s, double /*t*/, const double *x, const double *u, const double */*p*/) {

		const double c = 20;

		s(0, x_indexode(0)) = c*200.0 * x[0];
		s(0, x_indexode(1)) = c*200.0 * x[1];
		s(0, x_indexode(2)) = c*200.0 * x[2];
		s(0, x_indexode(3)) = c*200.0 * x[3];
		s(0, u_indexode(0)) = c*0.02 * u[0];
		s(0, u_indexode(1)) = c*0.02 * u[1];

		return true;
	}

	void ode(double *dx, double /*t*/, const double *x, const double *u, const double */*p*/) {

		double c = 20;
		dx[0] = (-10.0 * x[0] + u[0] + u[1]) * c;
		dx[1] = (-2.0 * x[1] + u[0] + 2.0 * u[1]) * c;
		dx[2] = (-3.0 * x[2] + 5.0 * x[3] + u[0] - u[1]) * c;
		dx[3] = (5.0 * x[2] - 3.0 * x[3] + u[0] + 3.0 * u[1]) * c;
	}

	bool ode_structure(tw::DiffStructure &s) {

		s(0, x_indexode(0));
		s(0, u_indexode(0));
		s(0, u_indexode(1));

		s(1, x_indexode(1));
		s(1, u_indexode(0));
		s(1, u_indexode(1));

		s(2, x_indexode(2));
		s(2, x_indexode(3));
		s(2, u_indexode(0));
		s(2, u_indexode(1));

		s(3, x_indexode(2));
		s(3, x_indexode(3));
		s(3, u_indexode(0));
		s(3, u_indexode(1));

		return true;
	}

	bool ode_diff(tw::DiffStructure& s, double /*t*/, const double* /*x*/, const double* /*u*/, const double* /*p*/) {

		const double c = 20.0;

		//dx[0] = (-10.0 * x[0] + u[0] + u[1]) * p[0];
		s(0, x_indexode(0)) = -10.0 * c;
		s(0, u_indexode(0)) = c;
		s(0, u_indexode(1)) = c;

		//dx[1] = (-2.0 * x[1] + u[0] + 2.0 * u[1]) * p[0];
		s(1, x_indexode(1)) = -2.0 * c;
		s(1, u_indexode(0)) = c;
		s(1, u_indexode(1)) = 2.0 * c;

		//dx[2] = (-3.0 * x[2] + 5.0 * x[3] + u[0] - u[1]) * p[0];
		s(2, x_indexode(2)) = -3.0 * c;
		s(2, x_indexode(3)) = 5.0 * c;
		s(2, u_indexode(0)) = c;
		s(2, u_indexode(1)) = -c;

		//dx[3] = (5.0 * x[2] - 3.0 * x[3] + u[0] + 3.0 * u[1]) * p[0];
		s(3, x_indexode(2)) = 5.0 * c;
		s(3, x_indexode(3)) = -3.0 * c;
		s(3, u_indexode(0)) = c;
		s(3, u_indexode(1)) = 3.0 * c;

		return true;
	}

	void var_boundary(double *x_low, double *x_upp) {

		x_low[x_index(0, 0)] = 2;
		x_upp[x_index(0, 0)] = 2;
		x_low[x_index(n_dis - 1, 0)] = 2;
		x_upp[x_index(n_dis - 1, 0)] = 2;

		x_low[x_index(0, 1)] = 1;
		x_upp[x_index(0, 1)] = 1;
		x_low[x_index(n_dis - 1, 1)] = 3;
		x_upp[x_index(n_dis - 1, 1)] = 3;

		x_low[x_index(0, 2)] = 2;
		x_upp[x_index(0, 2)] = 2;
		x_low[x_index(n_dis - 1, 2)] = 1;
		x_upp[x_index(n_dis - 1, 2)] = 1;

		x_low[x_index(0, 3)] = 1;
		x_upp[x_index(0, 3)] = 1;
		x_low[x_index(n_dis - 1, 3)] = -2;
		x_upp[x_index(n_dis - 1, 3)] = -2;

	}

	void neben(double *c, double t, const double *x, const double */*u*/, const double */*p*/) {
		const double cc = 20;

		c[0] = x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3] - 3.0*peak(t*cc,
			3, 12) - 3.0*peak(t*cc, 6, 10) - 3.0*peak(t*cc, 10, 6) -
			8.0*peak(t*cc, 15, 4);
		//c[0] = x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3] - 3.0*peak(t*p[0], 3, 12);
	}

	void neben_boundary(double *c_low, double */*c_upp*/) {
		c_low[0] = 0.01;
	}

	bool neben_structure(tw::DiffStructure &/*s*/) {

		//s(0, x_index(0, 0));
		//s(0, x_index(0, 1));
		//s(0, x_index(0, 2));
		//s(0, x_index(0, 3));
		return false;
	}

	bool neben_diff(tw::DiffStructure &/*s*/, double /*t*/, const double */*x*/, const double */*u*/, const double */*p*/) {

		// c[0] = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3] - 3.0*peak(t*p[0], 3, 12) - 3.0*peak(t*p[0], 6, 10) - 3.0*peak(t*p[0], 10, 6) - 8.0*peak(t*p[0], 15, 4);
		//s(0, x_index(0, 0)) = 2.0 * x[0];
		//s(0, x_index(0, 1)) = 2.0 * x[1];;
		//s(0, x_index(0, 2)) = 2.0 * x[2];;
		//s(0, x_index(0, 3)) = 2.0 * x[3];;

		return false;
	}
	
};

/////////////////////////////////////////////////////////////////////////////

int main(int argv, char* argc[]) {

	tw::TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv, argc);
	tw::TWfolder folder(&twparameter,0);

	AlpRiderPhase ph(twparameter.NDIS);
	ph.setSolver(&twparameter);
	folder.Add(&ph);
	
	tw::Viewer *viewer = nullptr;

	if (twparameter.PLOT) viewer = new tw::Viewer(&twparameter);

	folder.Init();
	folder.Init(viewer);
	
	folder.meshRef();

	delete viewer;

	return 0;
}
