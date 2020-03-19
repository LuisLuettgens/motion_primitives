/*----------------------------------------------------------------
 *
 * Brachistochrone 3 (2 Steuerungen)
 * M. Gerdts
 *
 *----------------------------------------------------------------*/

#include "TWfolder.h"
#include "TWproblem.h"
#include "TransWORHP.h"

double phasenPlot (int &len, int &ndgl, int &nsteuer, double *t,double *x, double *u, int &i, int &index);

class BrachistoPhase3 : public tw::TransWorhpProblem {
public:

	double grav;
	double c;
	
	BrachistoPhase3(int dis) : TransWorhpProblem(tw::TWdimension("Brachistochrone 3",dis,2,2,1,0,1)),
		grav(1),//9.81
		c(0.01) {
		freetime = true;
	}
	
	void OpenWindows(tw::Viewer *viewer) override {
		viewer->PhasePlot("Phasen-Plot", phasenPlot,0,1);
	}

	void x_init(double *x, int i, int dis) override {

		//x[0] = 0;
		//x[1] = 0;
	}

	void u_init(double *u, int i, int dis) override {

		u[0] = 1.0;
	}

	double obj() override {
		return p(0);
	}


	bool obj_structure(tw::DiffStructure &s) override {

		s(0,p_index(0));
		return true;
	}

	bool obj_diff_p(tw::DiffStructure &s) {

		s(0,p_index(0)) = 1.0;
		return true;
	}
	
	void ode(double *dx, double t, const double *x, const double *u, const double *p) override {

		dx[0] = sqrt( 2*grav*(c-x[1]) ) * u[0]*p[0];
		dx[1] = sqrt( 2*grav*(c-x[1]) ) * u[1]*p[0];
	}

	bool ode_structure(tw::DiffStructure &s) override {
		
		s(0,x_indexode(1));
		s(0,u_indexode(0));
		s(0,p_indexode(0));
		
		s(1,x_indexode(1));
		s(1,u_indexode(1));
		s(1,p_indexode(0));
		
		return true;
	}

	bool ode_diff(tw::DiffStructure &s, double t, const double *x, const double *u, const double *p) override {
		
		s(0,x_indexode(1)) = -grav*u[0]/sqrt(2*grav*(c-x[1]))*p[0];
		s(0,u_indexode(0)) = sqrt( 2*grav*(c-x[1]) ) *p[0];
		
		s(1,x_indexode(1)) = grav*u[0]/sqrt(2*grav*(c-x[1]))*p[0];
		s(1,u_indexode(1)) = sqrt( 2*grav*(c-x[1]) ) *p[0];
		
		return true;
	}
	
	bool ode_diff_p(tw::DiffStructure& s, double t, const double *x, const double *u, const double *p, int index) override {
		
		s(0,p_indexode(0)) = sqrt( 2*grav*(c-x[1]) ) * u[0];
		s(1,p_indexode(0)) = sqrt( 2*grav*(c-x[1]) ) * u[1];
		
		return true;
	}

	void u_boundary(double *u_low, double *u_upp) override {

		u_low[0] = -1.0;
		u_upp[0] =  1.0;
		u_low[1] = -1.0;
		u_upp[1] =  1.0;
	}

	void x_boundary(double *x_low, double *x_upp) override {

	}
	
	void p_boundary(double *p_low, double *p_upp) override {
		p_low[0] = 0.1;
		p_upp[0] = 10.0;
	}

	void var_boundary(double *x_low, double *x_upp) override {

		x_low[x_index(0,0)] = x_upp[x_index(0,0)] = 0.0;
		x_low[x_index(0,1)] = x_upp[x_index(0,1)] = 0.0;
		
		x_low[x_index(n_dis-1,0)] = x_upp[x_index(n_dis-1,0)] = 1.0;
		x_low[x_index(n_dis-1,1)] = x_upp[x_index(n_dis-1,1)] = -1.0;
	}
	
	void neben(double *c, double t, const double *x, const double *u, const double *p) override {

		c[0] = u[0]*u[0] + u[1]*u[1];
	}

	void neben_boundary(double *c_low, double *c_upp) override {
		c_low[0] = 1.0;
		c_upp[0] = 1.0;
	}
};

double phasenPlot (int &len, int &ndgl, int &nsteuer, double *t,double *x, double *u, int &i, int &index) {
	if (index==0) return x[i*(ndgl+nsteuer)+1];
	if (index==1) return x[i*(ndgl+nsteuer)];
}

/////////////////////////////////////////////////////////////////////////////

int main(int argv, char* argc[]) {

	tw::TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv,argc);
	tw::TWfolder folder(&twparameter,0);

	BrachistoPhase3 ph(twparameter.NDIS);
	ph.setSolver(&twparameter);
	folder.Add(&ph);

	tw::Viewer *viewer = nullptr;
	if (twparameter.PLOT) viewer = new tw::Viewer(&twparameter);

	folder.Init();
	folder.Init(viewer);

	folder.Loop(2);
	
	delete viewer;

	return 0;
}
