/*----------------------------------------------------------------
 *
 *  Tutorial: Splineproblem mit Mehrfachschiessverfahren
 *
 *----------------------------------------------------------------*/

#ifdef WIN32
#include "windows.h"
#endif

#include "ExplTransWORHP.h"


//#define RB

using namespace std;

class Spline_Expl : public ExplTransWorhp {
public:

	Spline_Expl(int dis, vector<int> &multinode) : 
#ifdef RB
ExplTransWorhp("Explicit Spline",dis,multinode, 3,1,0,2,0) {}
#else
ExplTransWorhp("Explicit Spline",dis,multinode, 3,1,0,0,0) {}
#endif

	double obj() {
		return x(n_dis-1,2);
	}

	bool obj_structure(DiffStructure &s) {
		s(0, x_index(n_dis-1,2));
		return true;
	}

	bool obj_diff(DiffStructure &s) {
		cout << "OBJ_" << endl;
		s(0, x_index(n_dis-1,2)) = 1;
		return true;
	}


	void ode(double *dx, double t, const double *x, const double *u,
			 const double *p) {

		dx[0] = x[1];
		dx[1] = u[0];
		dx[2] = u[0]*u[0];
	}

	bool ode_structure(DiffStructure &s) {
		s(0, x_indexode(1)); // dx[0] / dx[1]
		s(1, u_indexode(0)); // dx[1] / du[0]
		s(2, u_indexode(0)); // dx[2] / du[0]
		return true;
	}
/*
	bool ode_diff(DiffStructure &s, double t, const double *x,
			      const double *u, const double *p) {
return false;
		s(0, x_index(0,1))= 1;
		s(1, u_index(0,0))= 1;
		s(2, u_index(0,0))= 2*u[0];
		return true;
	}
*/
	void u_boundary(double *u_low, double *u_upp) {
		u_low[0] = -6;
		u_upp[0] = +6;
	}

	void x_boundary(double *x_low, double *x_upp) {
		x_low[2] = 0;
	}

	void var_boundary(double *x_low, double *x_upp) {
		x_low[x_index(0,0)] = x_upp[x_index(0,0)] = 0;
		x_low[x_index(0,1)] = x_upp[x_index(0,1)] = 1;
		x_low[x_index(0,2)] = x_upp[x_index(0,2)] = 0;

#ifdef RB
#else
		x_low[x_index(n_dis-1,0)] = x_upp[x_index(n_dis-1,0)] = 0;
		x_low[x_index(n_dis-1,1)] = x_upp[x_index(n_dis-1,1)] = 1;
#endif
	}

	void u_init(double *u, int i, int dis) {

		u[0] = -6 + (12.*i)/dis;
	}

	void x_init(double *x, int i, int dis) {

		x[0] = -2;
	}
	
#ifdef RB
	void rand(double *r) {

		r[0] = x(n_dis-1,0);
		r[1] = x(n_dis-1,1) - 1;
		
	}
#endif

/*	void neben_boundary(double *c_low, double *c_upp) {
		c_low[0] = -100; //-65./180.*M_PI;
		c_upp[0] = 100;
	}

	void neben(double *c, double t, const double *x, const double *u, const double *p) {
		c[0] = x[0]-x[1];		
	}

	bool neben_structure(DiffStructure &s) {
		s(0,x_index(0,0));
		s(0,x_index(0,1));
		return true;
	}
*/

};

///////////////////////////////////////////////////////////////////


// Aufruf: spline_expl
int main(int argv, char* argc[]) {

	TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv,argc);

	// Force explicit method... 
	twparameter.twdiscretization = TWdiscretization(TW_MultipleShooting,0,0);
	twparameter.USERHM=-1;
	//twparameter.USERDF=-1;
	twparameter.USERDG=-1;

	TWfolder folder(&twparameter,0);
	
	if (twparameter.multinode.size()==0) {
		twparameter.multinode.push_back(0);
		twparameter.multinode.push_back(twparameter.NDIS*.4);
		twparameter.multinode.push_back(twparameter.NDIS*.7);
		twparameter.multinode.push_back(twparameter.NDIS-1);
	}

	Spline_Expl ph(twparameter.NDIS, twparameter.multinode);
	ph.butcher->Init(twparameter.butchertableau,twparameter.stepsize);

	folder.Add(&ph);

	Viewer *viewer = 0;
	if (twparameter.PLOT) viewer = new Viewer(&twparameter);

	folder.Init();
	folder.Init(viewer);
	folder.Loop();

	delete viewer;

	return 0;
}

