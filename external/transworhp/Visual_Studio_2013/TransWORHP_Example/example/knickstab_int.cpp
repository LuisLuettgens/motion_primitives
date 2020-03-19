/*----------------------------------------------------------------
 *
 *  Example: Knickstab in Lagrange-Form
 *
 *----------------------------------------------------------------*/

#ifdef WIN32
#include "windows.h"
#endif

#include "TransWORHP.h"
#include "Viewer.h"
#include "TWfolder.h"
#include "base/twstatus.h"

using namespace std;

class KnickPhase : public TransWorhp {
public:
	
	double alpha;

	KnickPhase(int dis, double alpha) : TransWorhp("Knickstab", dis,2,1,0,0,0,1), alpha(alpha) {
	
		lagrange_weight[0] = 1;
	}

	string GetXTitle(int d) override {
		if (d==0) return "x";
		if (d==1) return "theta";
	}

	string GetUTitle(int d) override {
		if (d==0) return "u";
	}
	
	
/*	void x_init(double *x, int i, int dis) {
		x[0] = 0.05;
	}*/
	
	void u_init(double *u, int i, int dis) override {
		u[0] = 0.5;
	}
	

	void integral(double *f, double t, const double *x, const double *u,
			 const double *p) override {
		f[0] = 0.5 * u[0]*u[0] + alpha * cos(x[1]);
	}

	bool integral_structure(DiffStructure &s) override {
		s(0, u_indexode(0));
		s(0, x_indexode(1));
		return true;
	}
	
	bool integral_diff(DiffStructure &s, double t, const double *x, const double *u,
			const double *p) override {
	
		s(0, u_indexode(0)) = u[0];
		s(0, x_indexode(1)) = -alpha * sin(x[1]);
		return true;
	}
	
	double obj() override {
		return 0;
	}


	void ode(double *dx, double t, const double *x, const double *u, const double *p) override {
		dx[0] = sin(x[1]);
		dx[1] = u[0];
	}

	bool ode_structure(DiffStructure &s) override {
		s(0,x_index(0,1));
		s(1,u_index(0,0));
		return true;
	}
	
	bool ode_diff(DiffStructure &s, double t, const double *x, const double *u, const double *p) override {
		s(0,x_index(0,1)) = cos(x[1]);
		s(1,u_index(0,0)) = 1;
		return true;
	}
	
		
	void x_boundary(double *x_low, double *x_upp) override {
		x_low[0] = -.05;
		x_upp[0] = +.05;
	}

	void var_boundary(double *x_low, double *x_upp) override {
		x_low[x_index(0,0)] = x_upp[x_index(0,0)] = 0;
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

	Viewer *viewer = nullptr;
	if (twparameter.PLOT) viewer = new Viewer(&twparameter);

	folder.Init();
	folder.Init(viewer);
	
	folder.worhp_p.TolOpti = 1e-9;
	MyStatus("","Setze TolOpti auf " + ToString(folder.worhp_p.TolOpti),8 );
	
	folder.Loop();

	delete viewer;

	return 0;
}

