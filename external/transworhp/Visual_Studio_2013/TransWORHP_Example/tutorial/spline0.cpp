/*----------------------------------------------------------------
 *
 *  Tutorial: Splineproblem ohne Grafik
 *
 *----------------------------------------------------------------*/

#ifdef WIN32
#include "windows.h"
#endif

#include "TransWORHP.h"

using namespace std;

class Spline0 : public TransWorhp {
public:

	Spline0(int dis) : TransWorhp("Spline0",dis,3,1,0,0,0) {}

	double obj() {

		return x(n_dis-1,2);
	}

	void ode(double *dx, double t, const double *x, const double *u,
			 const double *p) {

		dx[0] = x[1];
		dx[1] = u[0];
		dx[2] = u[0]*u[0];
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

};

///////////////////////////////////////////////////////////////////

int main(int argv, char* argc[]) {

	TWparameter twparameter("transworhp.xml");
	TWfolder folder(&twparameter,0);

	Spline0 ph(11);
	folder.Add(&ph);

	folder.Init();

	folder.Loop();

	return 0;
}

