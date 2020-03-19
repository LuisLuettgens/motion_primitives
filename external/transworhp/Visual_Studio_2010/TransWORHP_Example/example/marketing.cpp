/*-----------------------------------------------------------------------
 *
 *
 *
 *-----------------------------------------------------------------------*/

#ifdef WIN32
#include "windows.h"
#endif

#include <iostream>
#include <iomanip>
#include <cmath>

#include "TransWORHP.h"

using namespace std;

class MarketingPhase : public TransWorhp {
public:

	MarketingPhase(int dis) : TransWorhp("Marketing",dis,3,1,0,0,0) {}

/*
	void init() {

		for (int i=0;i<n_dis;i++) {
			X[u_index(i,0)] = .5;
		}

	}*/
	
	void u_init(double *u, int i, int dis) {
	  
	u[0] = .5;
	
	  }

	double obj() {

		return  x(n_dis-1,2) * x(n_dis-1,2);
	}

	bool obj_structure(DiffStructure &s) {

		s(0, x_index(n_dis-1,2));
		return true;
	}
	
	bool obj_diff(DiffStructure &s) {
		s(0,x_index(n_dis-1,2)) = 2* x(n_dis-1,2);		
		return true;
	}
	

	void ode(double *dx, double t, const double *x, const double *u, const double *p) {

		dx[0] = x[1];
		dx[1] = u[0];
		dx[2] = u[0]*u[0];
	
	}

	// Optional
	bool ode_structure(DiffStructure &s) {

		s(0, x_indexode(1)); // dx[0] / dx[1]
		s(1, u_indexode(0)); // dx[1] / du[0]
		s(2, u_indexode(0)); // dx[2] / du[0]
		return true;
	}



	//void DG_diff_ode(double dg[100], double t, int var_index, int dis) {
	bool ode_diff(DiffStructure &s, double t, const double *x, const double *u, const double *p) {

		s(0, x_indexode(1))= 1;
		s(1, u_indexode(0))= 1;
		s(2, u_indexode(0))= 2*u[0];

		return true;
	}

	/*bool ode_diff_u(double *J, double t, const double *x, const double *u, const double *p, int index) {

	    if (index==0) { // nach u0
		J[0] = 0;
		J[1] = 1;
		J[2] = u[0];
	    }
	    return true;
	}*/

	bool ode_diff_p(DiffStructure &s, double t, const double *x, const double *u, const double *p, int index) {

		return true;
	}




	void u_boundary(double *u_low, double *u_upp) {

		u_low[0] = -6;
		u_upp[0] = +6;

	}

	void x_boundary(double *x_low, double *x_upp) {

		x_low[2] = 0;

	}

	void var_boundary(double *x_low, double *x_upp) {

		x_low[x_index ( 0,0 ) ] = x_upp[x_index ( 0,0 ) ] = 0;
		x_low[x_index ( 0,1 ) ] = x_upp[x_index ( 0,1 ) ] = 1;
		x_low[x_index ( 0,2 ) ] = x_upp[x_index ( 0,2 ) ] = 0;

		x_low[x_index ( n_dis-1,0 ) ] = x_upp[x_index ( n_dis-1,0 ) ] = 0;
		x_low[x_index ( n_dis-1,1 ) ] = x_upp[x_index ( n_dis-1,1 ) ] = 1;

	}

	/*
	 void rand(double *r) {

	  r[0] = x(0,0);//*x(0,0)*x(n_dis-1,0)*x(n_dis-1,0);
	  r[1] = x(0,1) - 1;
	  r[2] = x(0,2);

	  r[3] = x(n_dis-1,0);
	  r[4] = x(n_dis-1,1) - 1;

	  // r_boundary
	 }

	 bool rand_structure(DiffStructure &s) {

	  s(0,x_index(0,0)); //
	  s(0,x_index(n_dis-1,0)); //
	  s(1,x_index(0,1)); //
	  s(2,x_index(0,2)); //
	  s(3,x_index(n_dis-1,0)); //
	  s(4,x_index(n_dis-1,1)); //

	  return true;
	 }*/

};


/////////////////////////////////////////////////////////////////////////////


int main(int argv, char* argc[]) {

	XMLNode *xml = TransWorhp::ReadParams("transworhp.xml");

	int PLOT = 1;
	int NDIS = 11;

	for (int i=1;i<argv;i++) {
		if (strcmp(argc[i],"-n")==0) {
			i++;
			NDIS = atoi(argc[i]);
		}
		if (strcmp(argc[i],"+p")==0) {
			PLOT = 1;
		}

		if (strcmp(argc[i],"-p")==0) {
			PLOT = 0;
		}
	}

	/* Worhp data structures */
	OptVar    o;
	Workspace w;
	Params    p;
	Control   c;

	SplinePhase ph(NDIS);

	Viewer *viewer = 0;
	if (PLOT) viewer = new Viewer(xml);
	
	ph.Init(xml,o,w,p,c, viewer);
	ph.Loop();

	delete viewer;
	delete xml;
	
	return 0;
}
