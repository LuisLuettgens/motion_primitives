/*----------------------------------------------------------------
 *
 *  Example: Example 2 (Patterson, Hager, Rao, Benson, Huntington)
 *  unifedFrameworkAAS.pdf
 *
 *----------------------------------------------------------------*/

#ifdef WIN32
#include "windows.h"
#endif

#include "PmTransWORHP.h"

using namespace std;

class RaketenPhase : public PmTransWorhp {
  
private:
	double mu;
	double TT;
	double m0;
	double mDot;
  
public:
	
	RaketenPhase(TWdimension &TWdata) : PmTransWorhp(TWdata) {
		LinearTimeAxis(0,3.32);
		
		mu = 1;
		TT = 0.1405;
		m0 = 1;
		mDot = 0.0749;
	}

	void selectWindows(Viewer *viewer) override {

		viewer->AddStateView(0,"r");
		viewer->AddStateView(1,"theta");
		viewer->AddStateView(2,"v_r");
		viewer->AddStateView(3,"v_theta");
		
		viewer->AddControlView(0,"beta");
		
		viewer->AddMuView(0,"");
	}
	
	double obj() override {
		return -x(n_dis-1,0);
	}

	bool obj_structure(DiffStructure &s) override {
		s(0,x_index(n_dis-1,0));
		return true;
	}
	
	bool obj_diff(DiffStructure &s) override {
		s(0,x_index(n_dis-1,0)) = -1.;
		return true;
	}
	
	void x_init(double *x, int i, int dis) override {
		x[0] = 1.;
		x[1] = 1.;
		x[2] = 1.;
		x[3] = 1.;
	}
	
	void u_init(double *u, int i, int dis) override {
		u[0] = 1.;
	}
	
	void ode(double *dx, double t, const double *x, const double *u, const double *p) override {
		
		double a = TT/(m0 - mDot*t);
		
		dx[0] = x[2];
		dx[1] = x[3]/x[0];
		dx[2] = x[3]*x[3]/x[0] - mu/(x[0]*x[0]) + a*sin(u[0]);
		dx[3] = -x[2]*x[3]/x[0] + a*cos(u[0]);
	}

	bool ode_structure(DiffStructure &s) override {
		
		s(0,x_indexode(2));
		
		s(1,x_indexode(0));
		s(1,x_indexode(3));
		
		s(2,x_indexode(0));
		s(2,x_indexode(3));
		s(2,u_indexode(0));
		
		s(3,x_indexode(0));
		s(3,x_indexode(2));
		s(3,x_indexode(3));
		s(3,u_indexode(0));
		
		return true;
	}
	
	
	bool ode_diff(DiffStructure& s, double t, const double* x, const double* u, const double* p) override {
	  
		s(0,x_indexode(2)) = 1.;
		
		s(1,x_indexode(0)) = -x[3]/(x[0]*x[0]);
		s(1,x_indexode(3)) = 1./x[0];
		
		s(2,x_indexode(0)) = -(x[3]*x[3])/(x[0]*x[0]) - (-2*mu)/(x[0]*x[0]*x[0]);
		s(2,x_indexode(3)) = 2.*x[3];
		s(2,u_indexode(0)) = TT/(m0 - mDot*t) * cos(u[0]);
		
		s(3,x_indexode(0)) = -(-x[2]*x[3])/(x[0]*x[0]);
		s(3,x_indexode(2)) = -x[3]/x[0];
		s(3,x_indexode(3)) = -x[2]/x[0];
		s(3,u_indexode(0)) = TT/(m0 - mDot*t) * (-sin(u[0]));
		
		return true;
	}
	

	void rand(double *r) override {
	
		r[0] = x(n_dis-1,3) - sqrt(mu/x(n_dis-1,0));
	}
	
	bool rand_structure(DiffStructure &s) override {
	
		s(0,x_index(n_dis-1,0));
		s(0,x_index(n_dis-1,3));
	
		return true;
	}
	
	bool rand_diff(DiffStructure &s) override {
	
		s(0,x_index(n_dis-1,0)) = -(-sqrt(mu/x(n_dis-1,0)))/(2*x(n_dis-1,0));
		s(0,x_index(n_dis-1,3)) = 1.;
	
		return true;
	}

	void var_boundary(double *x_low, double *x_upp) override {
		x_low[x_index(0,0)] = x_upp[x_index(0,0)] = 1.0;
		x_low[x_index(0,1)] = x_upp[x_index(0,1)] = 0.0;
		x_low[x_index(0,2)] = x_upp[x_index(0,2)] = 0.0;
		x_low[x_index(0,3)] = x_upp[x_index(0,3)] = 1.0;
		
		x_low[x_index(n_dis-1,2)] = x_upp[x_index(n_dis-1,2)] = 0.0;
	}
	
};


/////////////////////////////////////////////////////////////////////////////

int main(int argv, char* argc[]) {

	TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv,argc);
	twparameter.twdiscretization = TWdiscretization(TW_Pseudospectral,0,0);
	TWfolder folder(&twparameter,0);
	
	TWdimension TWdim;
	TWdim.ID = "Orbit Raising";
	TWdim.n_dis = twparameter.NDIS;
	TWdim.n_ode = 4;
	TWdim.n_ctrl = 1;
	TWdim.n_param = 0;
	
	TWdim.n_neben = 0;
	TWdim.n_rand = 1;
	
	RaketenPhase ph(TWdim);
	folder.Add(&ph);

	Viewer *viewer = nullptr;
	if (twparameter.PLOT) viewer = new Viewer(&twparameter);
	
	folder.Init();
	folder.Init(viewer);
	
	folder.Loop();

	delete viewer;
	
	return 0;
}
