/*----------------------------------------------------------------
 *
 * Brachistochrone
 *
 *----------------------------------------------------------------*/

#ifdef WIN32
#include "windows.h"
#endif

#include "TransWORHP.h"

using namespace std;

class BrachistoPhase : public TransWorhp {
public:

	double grav;
	BrachistoPhase(int dis) : TransWorhp("Brachistochrone",dis,2,1,0,3,0) {
	
	grav = 9.81;
	}

	void x_init(double *x, int i, int dis) override {

		x[0] = -.1;

	}

	void u_init(double *u, int i, int dis) override {

		u[0] = -.2;

	}
	/*
	void init() {

	  for (int i=0;i<n_dis;i++) {
	  	X[u_index(i,0)] = -.2;
	  	X[x_index(i,0)] = -.1;
	  	//X[x_index(i,1)] =  + .3 * i/n_dis;
	  }

	}*/

	double obj() override {

		return x(n_dis-1,1);
	}


	bool obj_structure(DiffStructure &s) override {

		s(0,x_index(n_dis-1,1) );
		return true;

	}

	bool obj_diff(DiffStructure &s) override {

		s(0,x_index(n_dis-1,1)) = 1.;
		return true;

	}
	
	void ode(double *dx, double t, const double *x, const double *u, const double *p) override {

		dx[0] = u[0];
		dx[1] = sqrt( (1.+u[0]*u[0]) / (-2. * grav * x[0]) );
	}


	// Optional
	bool ode_structure(DiffStructure &s) override {
		//return false;
		
		s(0,u_index(0,0));

		s(1,u_index(0,0));
		s(1,x_index(0,0));

		return true;
	
	}

	bool ode_diff(DiffStructure &s, double t, const double *x, const double *u, const double *p) override {

		s(0,u_index(0,0)) = 1.;

		s(1,u_index(0,0)) = u[0]/sqrt(-(2+2*u[0]*u[0])*grav*x[0]) ;
		s(1,x_index(0,0)) = -0.5*(1+u[0]*u[0])/(sqrt(-(2+2*u[0]*u[0])*grav*x[0])*x[0]) ;
		
		return true;
		
	}

	void u_boundary(double *u_low, double *u_upp) override {

		// u_low[0] =  0;
		// u_upp[0] = 0;

	}

	void x_boundary(double *x_low, double *x_upp) override {

		//x_low[0] = 0;
		// x_upp[0] = 0;

	}

	void var_boundary(double *x_low, double *x_upp) override {

		/*r[0] = x(0,0) + 0.02;
		r[1] = x(0,1);
		r[2] = x(n_dis-1,0) + 1;
		*/

		/*x_low[x_index(0,0)] = -.02;
		x_upp[x_index(0,0)] = -.02;

		x_low[x_index(0,1)] = 0;
		x_upp[x_index(0,1)] = 0;*/

		//x_low[x_index(n_dis-1,0)] = -1;
		//x_upp[x_index(n_dis-1,0)] = -1;

	}

	void rand(double *r) override {
		r[0] = x(0,0) + 0.02;
		r[1] = x(0,1);
		r[2] = x(n_dis-1,0) + 1.;
	}


	bool rand_structure(DiffStructure &s) override {

		s(0,x_index(0,0));
		s(1,x_index(0,1));
		s(2,x_index(n_dis-1,0));
		return true;

	}
	
	bool rand_diff(DiffStructure &s) override {
	
		s(0,x_index(0,0)) = 1.; 
		s(1,x_index(0,1)) = 1.; 
		s(2,x_index(n_dis-1,0)) = 1.; 		
		return true;	
	}

};

/////////////////////////////////////////////////////////////////////////////

int main(int argv, char* argc[]) {

	TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv,argc);
	TWfolder folder(&twparameter,0);




	BrachistoPhase ph(twparameter.NDIS);
	folder.Add(&ph);

	Viewer *viewer = nullptr;
	if (twparameter.PLOT) viewer = new Viewer(&twparameter);

	folder.Init();
	folder.Init(viewer);
	
	folder.Loop();
	
	delete viewer;

	return 0;
}

