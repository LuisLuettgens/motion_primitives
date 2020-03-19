/*----------------------------------------------------------------
 *
 * Example: Chemischer Reaktor
 *
 *----------------------------------------------------------------*/

#ifdef WIN32
#include "windows.h"
#endif

#include "ExplTransWORHP.h"

using namespace std;

class ChemieExplPhase : public ExplTransWorhp {
public:

	ChemieExplPhase(int dis, vector<int> &multinode) : ExplTransWorhp("Chemiereaktor (expl)",dis,multinode,2,1,0,0,0) {}

	string GetXTitle(int d) override {
		if (d==0) return "x_1";
		if (d==1) return "x_2";
	}

	string GetUTitle(int d) override {
		if (d==0) return "u";
	}
	
	double obj() override {
		return -x(n_dis-1,1);
	}

	bool obj_structure(DiffStructure &s) override {
		//return false;
		s(0,x_index(n_dis-1,1) );
		return true;
	}
/*
	bool obj_diff(DiffStructure &s) {
		return false;
		s(0,x_index(n_dis-1,1) ) = -1;
		return true;
	}
*/	
	
	void ode(double *dx, double t, const double *x, const double *u, const double *p) override {
		dx[0] = -u[0]*x[0] +   u[0]*u[0]*x[1];
		dx[1] =  u[0]*x[0] - 3*u[0]*u[0]*x[1];
	}

	bool ode_structure(DiffStructure &s) override {
		s(0,u_indexode(0));
		s(0,x_indexode(0));
		s(0,x_indexode(1));

		s(1,u_indexode(0));
		s(1,x_indexode(0));
		s(1,x_indexode(1));
		return true;
	}
/*
	bool ode_diff(DiffStructure& s, double t, const double* x, const double* u, const double* p) {
		return false;
		s(0,u_index(0,0)) = -x[0] + 2*u[0]*x[1];
		s(0,x_index(0,0)) = -u[0];
		s(0,x_index(0,1)) = u[0]*u[0];

		s(1,u_index(0,0)) = x[0] - 6*u[0]*x[1];
		s(1,x_index(0,0)) = u[0];
		s(1,x_index(0,1)) = -3*u[0]*u[0];
		return true;
	}
*/
	
	bool ode_diff2(DiffStructure& s, double t, const double *x, const double *u, const double *p) override {
		s(0,0) = -u[0];
		s(0,1) = u[0]*u[0];
		s(1,0) = u[0];
		s(1,1) = -3*u[0]*u[0];
		return true;
	}

	void u_boundary(double *u_low, double *u_upp) override {
		u_low[0] =  0;
		u_upp[0] = +1;
	}

	void var_boundary(double *x_low, double *x_upp) override {
		x_upp[x_index(0,0)] = x_low[x_index(0,0)] = 1;
		x_upp[x_index(0,1)] = x_low[x_index(0,1)] = 0;
	}

};

/////////////////////////////////////////////////////////////////////////////

int main(int argv, char* argc[]) {

	TWparameter twparameter("transworhp.xml");
	std::map<std::string,std::string> A = twparameter.Arguments(argv,argc);
	int B =ToInt(A["div"]);
	if (B==0) B=3;
	
	cout << B << " " << twparameter.NDIS;
	if (B>twparameter.NDIS) return 0;
	
	// Force explicit method... 
	twparameter.twdiscretization = TWdiscretization(TW_MultipleShooting,0,0);
	twparameter.USERHM=-1;
	twparameter.USERDF=-1;
	
	TWfolder folder(&twparameter,0);
	
	vector<int> multinode;
	
	for (int i=0;i<B;i++) {
		int a = (double)(i * (twparameter.NDIS-1)) / (B-1);
		multinode.push_back(a);
	}
	
	ChemieExplPhase ph(twparameter.NDIS,multinode);
	ph.butcherInit(twparameter.butchertableau,twparameter.stepsize);
	
	folder.Add(&ph);
	
	
	Viewer *viewer = nullptr;
	if (twparameter.PLOT) viewer = new Viewer(&twparameter);
	
	folder.Init();
	folder.Init(viewer);
	
	folder.Loop(2);

	delete viewer;
	
	return 0;
}

