/*----------------------------------------------------------------
 *
 *  Example: Oszillator
 *
 *----------------------------------------------------------------*/

#ifdef WIN32
#include "windows.h"
#endif

#include "TransWORHP.h"

using namespace std;

class OszillatorPhase : public TransWorhp {
public:

	OszillatorPhase(int dis) : TransWorhp("Harmonischer Oszillator", dis,2,1,1,0,0) {
		freetime = 1;
	}


	string GetXTitle(int d) {
		if (d==0) return "x_1";
		if (d==1) return "x_2";
	}

	string GetUTitle(int d) {
		if (d==0) return "u";
	}


	void p_init(double *p) {
		p[0] = 5.2;
	}


	double obj() {
		return p(0);
	}

	bool obj_structure(DiffStructure &s) {
		s(0,p_index(0));
		return true;
	}

	bool obj_diff(DiffStructure &s) {
		s(0,p_index(0)) = 1;
		return true;
	}


	void ode(double *dx, double t, const double *x, const double *u, const double *p) {
		dx[0] = x[1] * p[0];
		dx[1] = (-x[0] + u[0]) * p[0];
	}

	bool ode_structure(DiffStructure &s) {
		s(0,x_indexode(1));
		s(0,p_indexode(0));

		s(1,x_indexode(0));
		s(1,u_indexode(0));
		s(1,p_indexode(0));
		return true;
	}

	bool ode_diff(DiffStructure &s, double t, const double *x, const double *u, const double *p) {
		s(0,x_indexode(1)) = p[0];
	
		s(1,x_indexode(0)) = -p[0];
		s(1,u_indexode(0)) = p[0];
		return true;
	}
	
	bool ode_diff_p(DiffStructure& s, double t, const double *x, const double *u, const double *p, int index) {
		s(0,p_indexode(0)) = x[1];
		s(1,p_indexode(0)) = -x[0] + u[0];
		return true;
	}


	void u_boundary(double *u_low, double *u_upp) {
		u_low[0] = -1;
		u_upp[0] = +1;
	}

	void p_boundary(double *p_low, double *p_upp) {
		p_low[0] = 1;
		p_upp[0] = 10;
	}

	void var_boundary(double *x_low, double *x_upp) {
		x_low[x_index(0,0)] = x_upp[x_index(0,0)] = 3;
		x_low[x_index(0,1)] = x_upp[x_index(0,1)] = 1;
		x_low[x_index(n_dis-1,0)] = x_upp[x_index(n_dis-1,0)] = 0;
		x_low[x_index(n_dis-1,1)] = x_upp[x_index(n_dis-1,1)] = 0;
	}

};



/////////////////////////////////////////////////////////////////////////////

int main(int argv, char* argc[]) {

	TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv,argc);
	TWfolder folder(&twparameter,0);

	OszillatorPhase ph(twparameter.NDIS);
	folder.Add(&ph);

	Viewer *viewer = 0;
	if (twparameter.PLOT) viewer = new Viewer(&twparameter);

	folder.Init();
	folder.Init(viewer);
	folder.Loop();

	delete viewer;

	return 0;
}

