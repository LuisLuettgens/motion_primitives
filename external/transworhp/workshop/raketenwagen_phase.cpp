/*----------------------------------------------------------------
 *
 *  Workshop: Raketenwagen in 2 Phasen
 *
 *----------------------------------------------------------------*/

#include "TWfolder.h"
#include "TWproblem.h"
#include "TransWORHP.h"

class RaketenPhase : public tw::TransWorhpProblem {
public:

	const int mode;

	RaketenPhase(const tw::TWdimension &TWdata, int mode_) : TransWorhpProblem(TWdata), mode(mode_) {}

	void p_init(double *p) override {
		p[0] = 1;
	}
	
	void x_init(double *x, int i, int dis) override {
		x[0] = 4;
		x[1] = -1;
	}

	double obj() override {
		return p(0);
	}
	
	bool obj_structure(tw::DiffStructure &s) override {
		s(0,p_index(0));

		return true;
	}

	bool obj_diff(tw::DiffStructure &s) override {
		s(0,p_index(0)) = 1;

		return true;
	}

	void ode(double *dx, double t, const double *x, const double *u, const double *p) override {

		double control = 0;

		if (mode == 0) control = -1;
		if (mode == 1) control = 1;

		dx[0] = x[1] * p[0];
		dx[1] = control * p[0];
	}

	// Optional
	bool ode_structure(tw::DiffStructure &s) override {
		s(0,x_indexode(1));
		s(0,p_indexode(0));

		s(1,p_indexode(0));

		return true;
	}


	bool ode_diff(tw::DiffStructure& s, double t, const double* x, const double* u, const double* p) override {
		s(0,x_indexode(1)) = p[0];

		return true;
	}
	
	bool ode_diff_p(tw::DiffStructure& s, double t, const double *x, const double *u, const double *p, int index) override {

		double control = 0;
		if (mode == 0) control = -1;
		if (mode == 1) control = 1;

		s(0,p_indexode(0)) = x[1];
		s(1,p_indexode(0)) = control;

		return true;
	}

	void p_boundary(double *p_low, double *p_upp) override {
		p_low[0] = .1;
		p_upp[0] = 10;
	}

	void var_boundary(double *x_low, double *x_upp) override {

		if (mode == 0) {
			x_low[x_index(0,0)] = x_upp[x_index(0,0)] = 4;
			x_low[x_index(0,1)] = x_upp[x_index(0,1)] = -1;
		}

		if (mode == 1) {
			x_low[x_index(n_dis-1,0)] = x_upp[x_index(n_dis-1,0)] = 0;
			x_low[x_index(n_dis-1,1)] = x_upp[x_index(n_dis-1,1)] = 0;
		}
	}
};

///////////////////////////////////////////////////////////////////

class RaketenFolder : public tw::TWfolder {
public:

	RaketenFolder(tw::TWparameter *p) : TWfolder(p,2) {}

	void g_boundary(double *g_low, double *g_upp) override {
			g_low[0] = 0;
			g_upp[0] = 0;
			
			g_low[1] = 0;
			g_upp[1] = 0;
	}

	void con(double *C) override {

		// zunaechst nur lineare Nebenbedingungen erlaubt!
		double *X = worhp_o.X;
		int index1 = phases[0]->x_index(phases[0]->n_dis-1,0) + phases[0]->solver->Delta1;
		int index2 = phases[1]->x_index(0,0)                  + phases[1]->solver->Delta1;

		C[0] = X[index1] - X[index2];
		C[1] = X[index1+1] - X[index2+1];
	}
	
	bool con_structure(tw::DiffStructure &s) override {

		int index1 = phases[0]->x_index(phases[0]->n_dis-1,0) + phases[0]->solver->Delta1;
		int index2 = phases[1]->x_index(0,0)                  + phases[1]->solver->Delta1;
		
		s(0,index1);
		s(0,index2);
		
		s(1,index1+1);
		s(1,index2+1);
		
		return true;
	}
	
	bool con_diff(tw::DiffStructure &s, int colindex) override {

		int index1 = phases[0]->x_index(phases[0]->n_dis-1,0) + phases[0]->solver->Delta1;
		int index2 = phases[1]->x_index(0,0)                  + phases[1]->solver->Delta1;
		
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

	tw::TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv,argc);
	
	tw::TWdimension TWdim1;
	TWdim1.ID = "Raketenwagen / Phase";
	TWdim1.n_dis = twparameter.NDIS;
	TWdim1.n_ode = 2;
	TWdim1.n_ctrl = 0; // Es wird der Umschaltpunkt bestimmt! deshalb n_ctrl=0
	TWdim1.n_param = 1;
	TWdim1.multinode = {0,TWdim1.n_dis/2,TWdim1.n_dis-1};
	
	tw::TWdimension TWdim2 = TWdim1;
	TWdim2.n_dis = TWdim1.n_dis/2 + 1;
	
	RaketenFolder folder(&twparameter);

	RaketenPhase ph1(TWdim1,0);
	RaketenPhase ph2(TWdim2,1);
	
	ph1.setSolver(&twparameter);
	ph2.setSolver(&twparameter);
	
	ph1.solver->LinearTimeAxis(0,1);
	ph2.solver->LinearTimeAxis(1,2);
	
	folder.Add(&ph1);
	folder.Add(&ph2);

	tw::Viewer *viewer = nullptr;
	if (twparameter.PLOT) viewer = new tw::Viewer(&twparameter);

	folder.Init();
	folder.Init(viewer);
	folder.Loop();

	delete viewer;

	return 0;
}
