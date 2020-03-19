/*----------------------------------------------------------------
 *
 *  Tutorial: Splineproblem mit zwei Phasen
 *
 *----------------------------------------------------------------*/

#ifdef WIN32
#include "windows.h"
#endif

#include "TransWORHP.h"

using namespace std;

class SplinePhase : public TransWorhp {
public:
	
	int mode;

	SplinePhase(int dis, int mode) : TransWorhp("Spline",dis,3,1,0,0,0), mode(mode) {}

	double obj() override {

		return x(n_dis-1,2);
	}

	bool obj_structure(DiffStructure &s) override {

		s(0, x_index(n_dis-1,2));
		return true;
	}

	bool obj_diff(DiffStructure &s) override {

		s(0, x_index(n_dis-1,2)) = 1;
		return true;
	}

	void ode(double *dx, double t, const double *x, const double *u,
			 const double *p) override {

		dx[0] = x[1];
		dx[1] = u[0];
		dx[2] = u[0]*u[0];
	}

	bool ode_structure(DiffStructure &s) override {

		s(0, x_indexode(1)); // dx[0] / dx[1]
		s(1, u_indexode(0)); // dx[1] / du[0]
		s(2, u_indexode(0)); // dx[2] / du[0]
		return true;
	}

	bool ode_diff(DiffStructure &s, double t, const double *x, const double *u, const double *p) override {

		s(0, x_indexode(1))= 1;
		s(1, u_indexode(0))= 1;
		s(2, u_indexode(0))= 2*u[0];
		return true;
	}

	void u_boundary(double *u_low, double *u_upp) override {

		u_low[0] = -6;
		u_upp[0] = +6;
	}

	void x_boundary(double *x_low, double *x_upp) override {

		x_low[2] = 0;
	}

	void var_boundary(double *x_low, double *x_upp) override {
		
		if (mode==0) x_low[x_index(0,0)] = x_upp[x_index(0,0)] = 0;
		if (mode==0) x_low[x_index(0,1)] = x_upp[x_index(0,1)] = 1;
		x_low[x_index(0,2)] = x_upp[x_index(0,2)] = 0;

		if (mode==1) x_low[x_index(n_dis-1,0)] = x_upp[x_index(n_dis-1,0)] = 0;
		if (mode==1) x_low[x_index(n_dis-1,1)] = x_upp[x_index(n_dis-1,1)] = 1;
	}

	void u_init(double *u, int i, int dis) override {

		u[0] = -6 + (12.*i)/dis;
	}

};


class SplineFolder : public TWfolder {
	
public:
	SplineFolder(TWparameter *p) : TWfolder(p,2) {}
	
	void g_boundary(double *x_low, double *x_upp) override {
	 
			x_low[0] = 0;
			x_upp[0] = 0;
			
			x_low[1] = 0;
			x_upp[1] = 0;

	 }
	 
	void con(double *C) override {
		
		// zunaechst nur lineare Nebenbedingungen erlaubt!
		double *X = worhp_o.X;
		int index1 = phases[0]->x_index(phases[0]->n_dis-1,0) + phases[0]->Delta1;
		int index2 = phases[1]->x_index(0,0)                  + phases[1]->Delta1;

		C[0] = X[index1] - X[index2];
		C[1] = X[index1+1] - X[index2+1];
	}
	
	bool con_structure(DiffStructure &s) override {

		int index1 = phases[0]->x_index(phases[0]->n_dis-1,0) + phases[0]->Delta1;
		int index2 = phases[1]->x_index(0,0)                  + phases[1]->Delta1;
		
		s(0,index1);
		s(0,index2);
		
		s(1,index1+1);
		s(1,index2+1);
		
		return true;
	}
	
	bool con_diff(DiffStructure &s, int colindex) override {

		int index1 = phases[0]->x_index(phases[0]->n_dis-1,0) + phases[0]->Delta1;
		int index2 = phases[1]->x_index(0,0)                  + phases[1]->Delta1;
		
		// C[0] = X[index1] - X[index2];
		s(0,index1) = +1;
		s(0,index2) = -1;
		
		//C[1] = X[index1+1] - X[index2+1];
		s(1,index1+1) = +1;
		s(1,index2+1) = -1;
		
		return true;
	}
	
};

///////////////////////////////////////////////////////////////////

int main(int argv, char* argc[]) {

	TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv, argc);
	SplineFolder folder(&twparameter);
	
	SplinePhase ph(twparameter.NDIS,0);
	SplinePhase ph2(twparameter.NDIS/2+1,1);
	
	ph.LinearTimeAxis(0,1/3.);
	ph2.LinearTimeAxis(1/3.,1);
	
	folder.Add(&ph);
	folder.Add(&ph2);

	Viewer *viewer = nullptr;
	if (twparameter.PLOT) viewer = new Viewer(&twparameter);

	folder.Init();
	folder.Init(viewer);
	folder.meshRef();

	delete viewer;

	return 0;
}

