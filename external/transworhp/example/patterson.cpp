/*----------------------------------------------------------------
 *
 *  Example: Example 1 (Patterson, Hager, Rao, Benson, Huntington)
 *  unifedFrameworkAAS.pdf
 * 
 * y*  = 4./(1+3.*exp(t))
 * u*  = y* /2
 * mu* = exp(2*ln(1+3*exp(t))-t)/(exp(-5)+6+9*exp(5))
 *
 *----------------------------------------------------------------*/

#include "TWfolder.h"
#include "TWproblem.h"
#include "TransWORHP.h"

class calcError : public tw::TransWorhpProblem {
public:
	
	calcError(const tw::TWdimension &TWdata) : TransWorhpProblem(TWdata) {
		
	}

	void localinit() override {
		solver->LinearTimeAxis(0,5);
	}

	void selectWindows(tw::Viewer *viewer) override {
		viewer->AddStateView(0,"y");
		viewer->AddControlView(0,"u");
		
		viewer->AddMuView(0,"y");
	}
	
	void x_init(double *x, int i, int dis) override {
		x[0] = 4./(1+3.*exp(solver->T[i]));
	}
	
	void u_init(double *u, int i, int dis) override {
		u[0] = (4./(1+3.*exp(solver->T[i])))/2.;
	}
	
	double obj() override {
		return -x(n_dis-1,0);
	}

	bool obj_structure(tw::DiffStructure &s) override {
		s(0,x_index(n_dis-1,0));
		return true;
	}
	
	bool obj_diff(tw::DiffStructure &s) override {
		s(0,x_index(n_dis-1,0)) = -1.;
		return true;
	}
	
	void ode(double *dx, double t, const double *x, const double *u, const double *p) override {
		dx[0] = -x[0]+x[0]*u[0]-u[0]*u[0];
	}

	bool ode_structure(tw::DiffStructure &s) override {
		s(0,x_indexode(0));
		s(0,u_indexode(0));
		return true;
	}

	bool ode_diff(tw::DiffStructure& s, double t, const double* x, const double* u, const double* p) override {
		s(0,x_indexode(0)) = -1. + u[0];
		s(0,u_indexode(0)) = x[0] - 2.*u[0];
		return true;
	}
	
	void var_boundary(double *x_low, double *x_upp) override {
		x_low[x_index(0,0)] = x_upp[x_index(0,0)] = 1.0;
	}
};


/////////////////////////////////////////////////////////////////////////////

int main(int argv, char* argc[]) {

	tw::TWparameter twparameter("transworhp.xml");
	auto args = twparameter.Arguments(argv, argc);

	tw::TWdimension TWdim;
	TWdim.ID = "Beispiel fuer spektrale Konvergenz";
	TWdim.n_ode = 1;
	TWdim.n_ctrl = 1;
	TWdim.n_dis = twparameter.NDIS;
	TWdim.setMultinodes(2);

	tw::TWfolder folder(&twparameter,0);

	calcError ph(TWdim);
	ph.setSolver(&twparameter);
	folder.Add(&ph);

	std::unique_ptr<tw::Viewer> viewer;
	if (twparameter.PLOT) viewer = std::unique_ptr<tw::Viewer>(new tw::Viewer(&twparameter));

	folder.Init();
	folder.Init(viewer.get());

	// Werte aus paper
	folder.worhp_p.TolOpti = 1e-15;
	folder.worhp_p.TolFeas = 2e-15;

	folder.Loop();

	return 0;
}
