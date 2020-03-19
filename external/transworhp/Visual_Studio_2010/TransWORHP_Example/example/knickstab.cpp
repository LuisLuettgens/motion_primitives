/*----------------------------------------------------------------
 *
 *  Example: Knickstab
 *
 *----------------------------------------------------------------*/

#ifdef WIN32
#include "windows.h"
#endif

#include "TransWORHP.h"
#include "Viewer.h"
#include "TWfolder.h"

using namespace std;

class KnickPhase : public TransWorhp {
public:
	
	double alpha;

	KnickPhase(int dis, double alpha) : TransWorhp("Knickstab", dis,3,1,0,0,0), alpha(alpha) {}

	string GetXTitle(int d) {
		if (d==0) return "x";
		if (d==1) return "theta";
		if (d==2) return "I";
	}

	string GetUTitle(int d) {
		if (d==0) return "u";
	}
	
	
	void x_init(double *x, int i, int dis) {
		x[0] = 0.05;
	}
	
	void u_init(double *u, int i, int dis) {
		u[0] = 0.5;
	}
	

	double obj() {
		return x(n_dis-1,2);
	}

	bool obj_structure(DiffStructure &s) {
		s(0,x_index(n_dis-1,2));
		return true;
	}

	bool obj_diff(DiffStructure &s) {
		s(0,x_index(n_dis-1,2)) = 1;
		return true;
	}
	

	void ode(double *dx, double t, const double *x, const double *u, const double *p) {
		dx[0] = sin(x[1]);
		dx[1] = u[0];
		dx[2] = 0.5 * u[0]*u[0] + alpha * cos(x[1]);
	}

	bool ode_structure(DiffStructure &s) {
		s(0,x_index(0,1));
		s(1,u_index(0,0));
		s(2,x_index(0,1));
		s(2,u_index(0,0));
		return true;
	}
	
	bool ode_diff(DiffStructure &s, double t, const double *x, const double *u, const double *p) {
		s(0,x_index(0,1)) = cos(x[1]);
		s(1,u_index(0,0)) = 1;
		s(2,x_index(0,1)) = -alpha * sin(x[1]);
		s(2,u_index(0,0)) = u[0];
		return true;
	}
	
		
	void x_boundary(double *x_low, double *x_upp) {
		x_low[0] = -.05;
		x_upp[0] = +.05;
	}

	void var_boundary(double *x_low, double *x_upp) {
		x_low[x_index(0,0)] = x_upp[x_index(0,0)] = 0;
		x_low[x_index(0,2)] = x_upp[x_index(0,2)] = 0;
		x_low[x_index(n_dis-1,0)] = x_upp[x_index(n_dis-1,0)] = 0;
	}

};



/////////////////////////////////////////////////////////////////////////////

int main(int argv, char* argc[]) {
	
	TWparameter twparameter("transworhp.xml");
	map<string,string> args = twparameter.Arguments(argv,argc);
	TWfolder folder(&twparameter,0);

	// Parameter alpha aus Kommandozeile lesen
	double alpha = (args["alpha"]!="") ? ToDouble(args["alpha"]) : 9.90027;
	
	KnickPhase ph(twparameter.NDIS, alpha);
	folder.Add(&ph);

	Viewer *viewer = 0;
	if (twparameter.PLOT) viewer = new Viewer(&twparameter);

	folder.Init();
	folder.Init(viewer);
	folder.Loop();

	delete viewer;

	return 0;
}

