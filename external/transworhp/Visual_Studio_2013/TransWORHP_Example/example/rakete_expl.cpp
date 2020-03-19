/*----------------------------------------------------------------
 *
 *  Example: Raketenwagen expl.
 *
 *----------------------------------------------------------------*/

#ifdef WIN32
#include "windows.h"
#endif

#include "ExplTransWORHP.h"

using namespace std;

class RaketenExplPhase : public ExplTransWorhp {
public:

	RaketenExplPhase(int dis, int multinode) : ExplTransWorhp("Raketenwagen (expl)", dis,multinode,2,1,1) {
		freetime = 1;
	}

	string GetXTitle(int d) override {
		if (d==0) return "x_1";
		if (d==1) return "x_2";
	}

	string GetUTitle(int d) override {
		if (d==0) return "u";
	}
	
	
	void p_init(double *p) override {
		p[0] = 3;
	}
	
	void u_init(double *u, int i, int dis) override {
		
		if (i < dis*0.5) {
			u[0] = -1;
		} else {
			u[0] = 1;
		}
	}
	
	double obj() override {
		return p(0);
	}

	bool obj_structure(DiffStructure &s) override {
		s(0,p_index(0) );
		return true;
	}
	
	bool obj_diff(DiffStructure &s) override {
		return false;
		s(0,p_index(0) ) = 1;
		return true;
	}
	

	void ode(double *dx, double t, const double *x, const double *u, const double *p) override {
		dx[0] = x[1] * p[0];
		dx[1] = u[0] * p[0];
	}

	bool ode_structure(DiffStructure &s) override {
		s(0,x_indexode(1));
		s(0,p_indexode(0));

		s(1,u_indexode(0));
		s(1,p_indexode(0));
		return true;
	}

	bool ode_diff(DiffStructure& s, double t, const double* x, const double* u, const double* p) override {
		return false;
		s(0,x_indexode(1)) = p[0];
		s(1,u_indexode(0)) = p[0];
		return true;
	}
	
	bool ode_diff_p(DiffStructure& s, double t, const double *x, const double *u, const double *p, int index) override {
		return false;
		s(0,p_indexode(0)) = x[1];
		s(1,p_indexode(0)) = u[0];
		return true;
	}
	
	
	bool ode_diff2(DiffStructure& s, double t, const double *x, const double *u, const double *p) override {
		s(0,x_indexode(1)) = p[0];
		return true;
	}
		
		
	void u_boundary(double *u_low, double *u_upp) override {
		u_low[0] = -1;
		u_upp[0] = +1;
	}

	void p_boundary(double *p_low, double *p_upp) override {
		p_low[0] = 1;
		p_upp[0] = 10;
	}

	void var_boundary(double *x_low, double *x_upp) override {
		x_low[x_index(0,0)] = x_upp[x_index(0,0)] = 4;
		x_low[x_index(0,1)] = x_upp[x_index(0,1)] = -1;
		x_low[x_index(n_dis-1,0)] = x_upp[x_index(n_dis-1,0)] = 0;
		x_low[x_index(n_dis-1,1)] = x_upp[x_index(n_dis-1,1)] = 0;
	}
	
	void neben(double *c, double t, const double *x, const double *u, const double *p) override {
		
	}
	
	void neben_boundary(double *c_low, double *c_upp) override {

	}

	bool neben_structure(DiffStructure &s) override {

		return false;
	}
	
	void terminate() override {
		
	}
	
};


/////////////////////////////////////////////////////////////////////////////

int main(int argv, char* argc[]) {

	TWparameter twparameter("transworhp.xml");
	std::map<std::string,std::string> A = twparameter.Arguments(argv,argc);
	int B = ToInt(A["div"]);
	
	// Force explicit method... 
	twparameter.twdiscretization = TWdiscretization(TW_MultipleShooting,0,0);
	twparameter.USERHM=-1;
	twparameter.USERDF=-1;
	//twparameter.USERDG=-1;
	
	TWfolder folder(&twparameter,0);
	
	RaketenExplPhase ph(twparameter.NDIS,B);
	
	ph.butcherInit(twparameter.butchertableau,twparameter.stepsize);
	folder.Add(&ph);
	
	Viewer *viewer = nullptr;
	if (twparameter.PLOT) viewer = new Viewer(&twparameter);
	
	folder.Init();
	folder.Init(viewer);
	
	folder.Loop();
	
	delete viewer;
	
	return 0;
}

