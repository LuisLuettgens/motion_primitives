/*----------------------------------------------------------------
 *
 *  Example: Knickstab
 *
 *----------------------------------------------------------------*/

#include "TWfolder.h"
#include "TWproblem.h"
#include "TransWORHP.h"

class KnickPhase : public tw::TransWorhpProblem {
public:
	
	const double alpha;

	KnickPhase(int dis, double alpha) : TransWorhpProblem(tw::TWdimension("Knickstab", dis,3,1,0,0,0)), alpha(alpha) {}

	void selectWindows(tw::Viewer *viewer) override {

		viewer->AddStateView(0,"x");
		viewer->AddStateView(1,"theta");
		viewer->AddStateView(2,"I");
		
		viewer->AddControlView(0,"u");
	}
	
	void x_init(double *x, int i, int dis) override {
		x[0] = 0.05;
	}
	
	void u_init(double *u, int i, int dis) override {
		u[0] = 0.5;
	}
	
	double obj() override {
		return x(n_dis-1,2);
	}

	bool obj_structure(tw::DiffStructure &s) override {
		s(0,x_index(n_dis-1,2));
		return true;
	}

	bool obj_diff(tw::DiffStructure &s) override {
		s(0,x_index(n_dis-1,2)) = 1;
		return true;
	}
	
	void ode(double *dx, double t, const double *x, const double *u, const double *p) override {
		dx[0] = sin(x[1]);
		dx[1] = u[0];
		dx[2] = 0.5 * u[0]*u[0] + alpha * cos(x[1]);
	}

	bool ode_structure(tw::DiffStructure &s) override {
		s(0,x_index(0,1));
		s(1,u_index(0,0));
		s(2,x_index(0,1));
		s(2,u_index(0,0));
		return true;
	}
	
	bool ode_diff(tw::DiffStructure &s, double t, const double *x, const double *u, const double *p) override {
		s(0,x_index(0,1)) = cos(x[1]);
		s(1,u_index(0,0)) = 1;
		s(2,x_index(0,1)) = -alpha * sin(x[1]);
		s(2,u_index(0,0)) = u[0];
		return true;
	}
	
	void x_boundary(double *x_low, double *x_upp) override {
		x_low[0] = -.05;
		x_upp[0] = +.05;
	}

	void var_boundary(double *x_low, double *x_upp) override {
		x_low[x_index(0,0)] = x_upp[x_index(0,0)] = 0.0;
		x_low[x_index(0,2)] = x_upp[x_index(0,2)] = 0.0;
		x_low[x_index(n_dis-1,0)] = x_upp[x_index(n_dis-1,0)] = 0.0;
	}
};

/////////////////////////////////////////////////////////////////////////////

int main(int argv, char* argc[]) {

	tw::TWparameter twparameter("transworhp.xml");
	std::map<std::string,std::string> args = twparameter.Arguments(argv,argc);
	tw::TWfolder folder(&twparameter,0);

	// Parameter alpha aus Kommandozeile lesen
	double alpha = (args["alpha"]!="") ? std::stod(args["alpha"]) : 9.90027;

	KnickPhase ph(twparameter.NDIS, alpha);
	ph.setSolver(&twparameter);
	folder.Add(&ph);

	tw::Viewer *viewer = nullptr;
	if (twparameter.PLOT) viewer = new tw::Viewer(&twparameter);

	folder.Init();
	folder.Init(viewer);
	folder.Loop();

	delete viewer;

	return 0;
}
