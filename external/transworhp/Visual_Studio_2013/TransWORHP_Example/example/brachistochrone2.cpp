/*----------------------------------------------------------------
 *
 * Brachistochrone 2
 * M. Gerdts
 *
 *----------------------------------------------------------------*/

#ifdef WIN32
#include "windows.h"
#endif

#include "TransWORHP.h"
#include "../src/core/TWparameter.h"
#include "../src/core/Viewer.h"

using namespace std;

double phasenPlot (int &len, int &ndgl, int &nsteuer, double *t,double *x, double *u, int &i, int &index);

class BrachistoPhase2 : public TransWorhp {
public:

	double grav;
	double c;
	BrachistoPhase2(int dis) : TransWorhp("Brachistochrone 2",dis,2,1,1,0,0) {
	
	freetime = 1;
	  
	grav = 1;//9.81;
	c =0.01;
	
	}
	
	void OpenWindows(Viewer *viewer) override {
		viewer->PhasePlot("Phasen-Plot", phasenPlot,0,1);
	}

	void x_init(double *x, int i, int dis) override {

		//x[0] = 0;
		//x[1] = 0;
	}

	void u_init(double *u, int i, int dis) override {

		u[0] = 1;
	}

	double obj() override {
		return p(0);
	}


	bool obj_structure(DiffStructure &s) override {

		s(0,p_index(0));
		return true;
	}

	bool obj_diff_p(DiffStructure &s) {

		s(0,p_index(0)) = 1;
		return true;
	}
	
	void ode(double *dx, double t, const double *x, const double *u, const double *p) override {

		dx[0] = sqrt( 2*grav*(c-x[1]) ) * cos(u[0])*p[0];
		dx[1] = sqrt( 2*grav*(c-x[1]) ) * sin(u[0])*p[0];
	}

	bool ode_structure(DiffStructure &s) override {
		
		s(0,x_indexode(1));
		s(0,u_indexode(0));
		s(0,p_indexode(0));
		
		s(1,x_indexode(1));
		s(1,u_indexode(0));
		s(1,p_indexode(0));
		
		return true;
	}

	bool ode_diff(DiffStructure &s, double t, const double *x, const double *u, const double *p) override {
		
		s(0,x_indexode(1)) = -grav*cos(u[0])/sqrt(2*grav*(c-x[1]))*p[0];
		s(0,u_indexode(0)) = -sqrt(2*grav*(c-x[1]))*sin(u[0])*p[0];
		
		s(1,x_indexode(1)) = grav*sin(u[0])/sqrt(2*grav*(c-x[1]))*p[0];
		s(1,u_indexode(0)) = sqrt(2*grav*(c-x[1]))*cos(u[0])*p[0];
		
		return true;
	}
	
	bool ode_diff_p(DiffStructure& s, double t, const double *x, const double *u, const double *p, int index) override {
	
		s(0,p_indexode(0)) = sqrt( 2*grav*(c-x[1]) ) * cos(u[0]);
		s(1,p_indexode(0)) = sqrt( 2*grav*(c-x[1]) ) * sin(u[0]);
		
		return true;
	}

	void u_boundary(double *u_low, double *u_upp) override {

		//u_low[0] = -2*M_PI;
		//u_upp[0] =  2*M_PI;
	}

	void x_boundary(double *x_low, double *x_upp) override {

		x_upp[1] = 0;

	}
	
	void p_boundary(double *p_low, double *p_upp) override {
		p_low[0] = 0.1;
		p_upp[0] = 10;
	}

	void var_boundary(double *x_low, double *x_upp) override {

		x_low[x_index(0,0)] = x_upp[x_index(0,0)] = 0;
		x_low[x_index(0,1)] = x_upp[x_index(0,1)] = 0;
		
		x_low[x_index(n_dis-1,0)] = x_upp[x_index(n_dis-1,0)] = 1;
		x_low[x_index(n_dis-1,1)] = x_upp[x_index(n_dis-1,1)] = -1;
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




	BrachistoPhase2 ph(twparameter.NDIS);
	folder.Add(&ph);

	Viewer *viewer = nullptr;
	if (twparameter.PLOT) viewer = new Viewer(&twparameter);

	folder.Init();
	folder.Init(viewer);

	
	folder.Loop(2);
	
	delete viewer;

	return 0;
}

