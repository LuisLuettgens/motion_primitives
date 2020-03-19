/*----------------------------------------------------------------
 *
 * Example: Chemischer Reaktor
 *
 *----------------------------------------------------------------*/

#ifdef WIN32
#include "windows.h"
#endif

#include "TransWORHP.h"

using namespace std;

class ChemiePhase : public TransWorhp {
public:

	ChemiePhase(int dis) : TransWorhp("Chemiereaktor",dis,2,1,0,0,0) {}

	string GetXTitle(int d) {
		if (d==0) return "x_1";
		if (d==1) return "x_2";
	}

	string GetUTitle(int d) {
		if (d==0) return "u";
	}
	
	
	double obj() {
		return -x(n_dis-1,1);
	}

	bool obj_structure(DiffStructure &s) {
		s(0,x_index(n_dis-1,1) );
		return true;
	}

	bool obj_diff(DiffStructure &s) {
		s(0,x_index(n_dis-1,1) ) = -1;
		return true;
	}
	
	
	void ode(double *dx, double t, const double *x, const double *u, const double *p) {
		dx[0] = -u[0]*x[0] +   u[0]*u[0]*x[1];
		dx[1] =  u[0]*x[0] - 3*u[0]*u[0]*x[1];
	}

	bool ode_structure(DiffStructure &s) {
		s(0,u_index(0,0));
		s(0,x_index(0,0));
		s(0,x_index(0,1));

		s(1,u_index(0,0));
		s(1,x_index(0,0));
		s(1,x_index(0,1));
		return true;
	}

	bool ode_diff(DiffStructure& s, double t, const double* x, const double* u, const double* p) {
		s(0,u_index(0,0)) = -x[0] + 2*u[0]*x[1];
		s(0,x_index(0,0)) = -u[0];
		s(0,x_index(0,1)) = u[0]*u[0];

		s(1,u_index(0,0)) = x[0] - 6*u[0]*x[1];
		s(1,x_index(0,0)) = u[0];
		s(1,x_index(0,1)) = -3*u[0]*u[0];
		return true;
	}


	void u_boundary(double *u_low, double *u_upp) {
		u_low[0] =  0;
		u_upp[0] = +1;
	}

	void var_boundary(double *x_low, double *x_upp) {
		x_upp[x_index(0,0)] = x_low[x_index(0,0)] = 1;
		x_upp[x_index(0,1)] = x_low[x_index(0,1)] = 0;
	}

};

/////////////////////////////////////////////////////////////////////////////

int main(int argv, char* argc[]) {

	TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv,argc);
	TWfolder folder(&twparameter,0);

	ChemiePhase ph(twparameter.NDIS);
	folder.Add(&ph);

	Viewer *viewer = 0;
	if (twparameter.PLOT) viewer = new Viewer(&twparameter);

	folder.Init();
	folder.Init(viewer);
	folder.Loop();

	delete viewer;

	return 0;
}

