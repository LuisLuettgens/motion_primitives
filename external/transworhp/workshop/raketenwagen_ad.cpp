/*----------------------------------------------------------------
 *
 *  Tutorial: Raketenwagen
 *
 *----------------------------------------------------------------*/

#ifdef WIN32
#include "windows.h"
#endif

#include "MagicTransWORHP.h"

using namespace std;

class RaketenAD : public MagicTransWorhp {
public:

	RaketenAD(int dis) : MagicTransWorhp("Raketenwagen / Automatische Differentiation",dis,2,1,1,0,0) {
	
		freetime = 1;
	}

	void p_init(double *p) {

		p[0] = 3;

	}


	void obj(MagicDouble *F, MagicDouble *X) {
	
		F[0] = X[p_index(0)];
	}
 
 	void ode(MagicDouble *dx, double t, const MagicDouble *x, const MagicDouble *u,
			 const MagicDouble *p) {

		dx[0] = x[1] * p[0];
		dx[1] = u[0] * p[0];
	}

	
	void u_boundary(double *u_low, double *u_upp) {

		u_low[0] = -1;
		u_upp[0] = +1;

	}

	void x_boundary(double *x_low, double *x_upp) {

	}

	void p_boundary(double *p_low, double *p_upp) {
		p_low[0] = 1;
		p_upp[0] = 10;
	}

	void var_boundary(double *x_low, double *x_upp) {

		x_low[x_index(0,0)] = x_upp[x_index(0,0)] = 4;
		x_low[x_index(0,1)] = x_upp[x_index(0,1)] = -1;
		x_low[x_index(n_dis-1,0)] = x_upp[x_index(n_dis-1,0)] = 0;
		x_low[x_index(n_dis-1,1)] = x_upp[x_index(n_dis-1,1)] = 0;
		
	}


};


///////////////////////////////////////////////////////////////////

int main(int argv, char* argc[]) {

	TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv,argc);
	TWfolder folder(&twparameter,0);

	RaketenAD ph(twparameter.NDIS);
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

