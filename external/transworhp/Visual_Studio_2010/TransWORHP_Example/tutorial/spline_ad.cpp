/*----------------------------------------------------------------
 *
 *  Tutorial: Splineproblem mit automatischer Differentiation
 *
 *----------------------------------------------------------------*/

#ifdef WIN32
#include "windows.h"
#endif

#include "MagicTransWORHP.h"

using namespace std;

class SplineAD : public MagicTransWorhp {
public:

	SplineAD(int dis) : MagicTransWorhp("Spline / Automatische Differentiation",dis,3,1,0,0,0) {}

	void obj(MagicDouble *F, MagicDouble *X) {
		
		F[0] = X[x_index(n_dis-1,2)];
	}
 
	void ode(MagicDouble *dx, double t, const MagicDouble *x, const MagicDouble *u,
			 const MagicDouble *p) {

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

	void u_init(double *u, int i, int dis) {

		u[0] = -6 + (12.*i)/dis;
	}

};

///////////////////////////////////////////////////////////////////

int main(int argv, char* argc[]) {

	TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv,argc);
	TWfolder folder(&twparameter,0);

	SplineAD ph(twparameter.NDIS);
	folder.Add(&ph);

	Viewer *viewer = 0;
	if (twparameter.PLOT) viewer = new Viewer(&twparameter);

	folder.Init();
	folder.Init(viewer);
	folder.Loop();

	MagicDouble::Info();
	delete viewer;

	return 0;
}

