/*----------------------------------------------------------------
 *
 *  Example: Splineproblem
 *
 *----------------------------------------------------------------*/

#ifdef WIN32
#include "windows.h"
#endif

#include "TransWORHP.h"

using namespace std;

class SplinePhase : public TransWorhp {
public:

	double alpha;
	
	SplinePhase(int dis, double alpha) : TransWorhp("Spline",dis,3,1,0,0,0), alpha(alpha) {}

	string GetXTitle(int d) override {
		if (d==0) return "x_1";
		if (d==1) return "x_2";
		if (d==2) return "x_3";
	}

	string GetUTitle(int d) override {
		if (d==0) return "u";
	}
	

	double obj() override {
		return .5 * x(n_dis-1,2);
	}

	bool obj_structure(DiffStructure &s) override {
		s(0, x_index(n_dis-1,2));
		return true;
	}

	bool obj_diff(DiffStructure &s) override {
		s(0,x_index(n_dis-1,2)) = .5;
		return true;
	}


	void ode(double *dx, double t, const double *x, const double *u, const double *p) override {
		dx[0] = x[1];
		dx[1] = u[0];
		dx[2] = u[0]*u[0];
	}

	bool ode_structure(DiffStructure &s) override {
		s(0, x_indexode(1)); // dx[0] / dx[1]
		s(1, u_indexode(0)); // dx[1] / du[0]
		s(2, u_indexode(0)); // dx[2] / du[0]
		return true;
	}

	bool ode_diff(DiffStructure &s, double t, const double *x, const double *u, const double *p) override {
		s(0, x_indexode(1))= 1;
		s(1, u_indexode(0))= 1;
		s(2, u_indexode(0))= 2*u[0];
		return true;
	}


	void x_boundary(double *x_low, double *x_upp) override {
		x_upp[0] = alpha;
	}

	void var_boundary(double *x_low, double *x_upp) override {
		x_low[x_index ( 0,0 ) ] = x_upp[x_index ( 0,0 ) ] = 0;
		x_low[x_index ( 0,1 ) ] = x_upp[x_index ( 0,1 ) ] = 1;
		x_low[x_index ( 0,2 ) ] = x_upp[x_index ( 0,2 ) ] = 0;

		x_low[x_index ( n_dis-1,0 ) ] = x_upp[x_index ( n_dis-1,0 ) ] = 0;
		x_low[x_index ( n_dis-1,1 ) ] = x_upp[x_index ( n_dis-1,1 ) ] = -1;
	}

};

/////////////////////////////////////////////////////////////////////////////

int main(int argv, char* argc[]) {

	TWparameter twparameter("transworhp.xml");
	map<string,string> args = twparameter.Arguments(argv,argc);
	TWfolder folder(&twparameter,0);

	// Parameter alpha aus Kommandozeile lesen
	double alpha = (args["alpha"]!="") ? ToDouble(args["alpha"]) : 1/6.;
	
	SplinePhase ph(twparameter.NDIS, alpha);
	folder.Add(&ph);

	Viewer *viewer = nullptr;
	if (twparameter.PLOT) viewer = new Viewer(&twparameter);

	folder.Init();
	folder.Init(viewer);
	folder.Loop();

	delete viewer;

	return 0;
}

