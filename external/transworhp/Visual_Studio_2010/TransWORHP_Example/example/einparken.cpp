/*----------------------------------------------------------------
 *
 *  Example: Einparken
 *
 *----------------------------------------------------------------*/

#ifdef WIN32
#include "windows.h"
#endif

#include "TransWORHP.h"

using namespace std;

double phasenPlot (int &len, int &ndgl, int &nsteuer, double *t,double *x, double *u, int &i, int &index);

class EinparkenPhase : public TransWorhp {
public:

	EinparkenPhase(int dis) : TransWorhp("Einparken", dis,6,2,1,0,0) {
		freetime = 1;
	}


	void OpenWindows(Viewer *viewer) {
		viewer->PhasePlot("Phasen-Plot", phasenPlot,0,1);
	}

	string GetXTitle(int d) {
		if (d==0) return "x";
		if (d==1) return "y";
		if (d==2) return "v";
		if (d==3) return "theta";
		if (d==4) return "phi";
		if (d==5) return "hilfe";
	}

	string GetUTitle(int d) {
		if (d==0) return "Beschleunigung";
		if (d==1) return "sigma";
	}
	
	
	void p_init(double *p) {
		p[0] = 2;
	}

	
	double obj() {
		return p(0) + 0.0*x(n_dis-1, 5);
	}

	bool obj_structure(DiffStructure &s) {
		s(0,x_index(n_dis-1,5));
		s(0,p_index(0));
		return true;
	}
	
	bool obj_diff(DiffStructure &s) {
		s(0,x_index(n_dis-1,5)) = 0.0;
		s(0,p_index(0)) = 1;
		return true;
	}
	

	void ode(double *dx, double t, const double *x, const double *u, const double *p) {
		dx[0] = x[2] * cos(x[3]) * p[0];
		dx[1] = x[2] * sin(x[3]) * p[0];
		dx[2] = u[0] * p[0];
		dx[3] = x[2] * sin(x[4]) * p[0] / 1; // b = 1
		dx[4] = u[1] * p[0];
		dx[5] = u[0]*u[0] + u[1]*u[1];
	}

	bool ode_structure(DiffStructure &s) {
		s(0,x_indexode(2));
		s(0,x_indexode(3));
		s(0,p_indexode(0));

		s(1,x_indexode(2));
		s(1,x_indexode(3));
		s(1,p_indexode(0));

		s(2,u_indexode(0));
		s(2,p_indexode(0));

		s(3,x_indexode(2));
		s(3,x_indexode(4));
		s(3,p_indexode(0));

		s(4,u_indexode(1));
		s(4,p_indexode(0));

		s(5,u_indexode(0));
		s(5,u_indexode(1));
		return true;
	}

	bool ode_diff(DiffStructure& s, double t, const double* x, const double* u, const double* p) {
		s(0,x_indexode(2)) = cos(x[3]) * p[0];
		s(0,x_indexode(3)) = -x[2] * sin(x[3]) * p[0];

		s(1,x_indexode(2)) = sin(x[3]) * p[0];
		s(1,x_indexode(3)) = x[2] * cos(x[3]) * p[0];		

		s(2,u_indexode(0)) = p[0];

		s(3,x_indexode(2)) = sin(x[4]) * p[0] / 1; // b = 1
		s(3,x_indexode(4)) = x[2] * cos(x[4]) * p[0] / 1; // b = 1

		s(4,u_indexode(1)) = p[0];

		s(5,u_indexode(0)) = 2*u[0];
		s(5,u_indexode(1)) = 2*u[1];
		return true;
	}
	
	bool ode_diff_p(DiffStructure& s, double t, const double *x, const double *u, const double *p, int index) {
		s(0,p_indexode(0)) = x[2] * cos(x[3]);
		s(1,p_indexode(0)) = x[2] * sin(x[3]);
		s(2,p_indexode(0)) = u[0];
		s(3,p_indexode(0)) = x[2] * sin(x[4]) / 1; // b = 1
		s(4,p_indexode(0)) = u[1];
		return true;
	}
	
	void x_boundary(double *x_low, double *x_upp) {
		x_low[2] = 0;
		x_upp[2] = 10;
		x_low[4] = -M_PI/2; // phi_max
		x_upp[4] = M_PI/2;
		x_low[5] = 0;
	}
	
	void u_boundary(double *u_low, double *u_upp) {
		u_low[0] = -10;
		u_upp[0] = 10;
		u_low[1] = -M_PI/2; // sigma_max
		u_upp[1] = M_PI/2;
	}

	void p_boundary(double *p_low, double *p_upp) {
		p_low[0] = 1;
		p_upp[0] = 10;
	}

	void var_boundary(double *x_low, double *x_upp) {
		x_low[x_index(0,0)] = x_upp[x_index(0,0)] = 0;
		x_low[x_index(0,1)] = x_upp[x_index(0,1)] = 0;
		x_low[x_index(0,2)] = x_upp[x_index(0,2)] = 0;
		x_low[x_index(0,3)] = x_upp[x_index(0,3)] = 0;
		x_low[x_index(0,4)] = x_upp[x_index(0,4)] = 0;
		x_low[x_index(n_dis-1,0)] = x_upp[x_index(n_dis-1,0)] = 5;
		x_low[x_index(n_dis-1,1)] = x_upp[x_index(n_dis-1,1)] = -2;
		x_low[x_index(n_dis-1,2)] = x_upp[x_index(n_dis-1,2)] = 0;
		x_low[x_index(n_dis-1,3)] = x_upp[x_index(n_dis-1,3)] = -M_PI/2;
		x_low[x_index(n_dis-1,4)] = x_upp[x_index(n_dis-1,4)] = 0;
	}

};

double phasenPlot (int &len, int &ndgl, int &nsteuer, double *t,double *x, double *u, int &i, int &index) {
	if (index==0) return x[i*(ndgl+nsteuer)+1];
	if (index==1) return x[i*(ndgl+nsteuer)];
}

/////////////////////////////////////////////////////////////////////////////

int main(int argv, char* argc[]) {

	TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv,argc);
	TWfolder folder(&twparameter,0);

	EinparkenPhase ph(twparameter.NDIS);
	folder.Add(&ph);

	Viewer *viewer = 0;
	if (twparameter.PLOT) viewer = new Viewer(&twparameter);
	
	folder.Init();
	folder.Init(viewer);
	folder.Loop();

	delete viewer;
	
	return 0;
}

