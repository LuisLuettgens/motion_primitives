/*----------------------------------------------------------------
 *
 * Brachistochrone
 *
 * Beispiel 1
 * aus "Convergence of a Gauss Pseudospectral Method for Optimal Control"
 * von Hou, Hager, Rao
 *
 *----------------------------------------------------------------*/

#include "TWfolder.h"
#include "TWproblem.h"
#include "TransWORHP.h"

class BrachistoPhase : public tw::TransWorhpProblem {
public:

	const double grav;
	
	BrachistoPhase(const tw::TWdimension &TWdim) : TransWorhpProblem(TWdim), grav(10) {
		freetime = true;
	}

	void x_init(double *x, int i, int dis) override {
		
	}

	void u_init(double *u, int i, int dis) override {
		
	}

	void p_init(double *p) override {
		
	}

	double obj() override {
		return p(0);
	}

	bool obj_structure(tw::DiffStructure &s) override {
		s(0,p_index(0));
		return true;
	}

	bool obj_diff(tw::DiffStructure &s) override {
		s(0,p_index(0)) = 1.;
		return true;
	}
	
	void ode(double *dx, double t, const double *x, const double *u, const double *p) override {
		dx[0] = x[2]*sin(u[0]);
		dx[1] = -x[2]*cos(u[0]);
		dx[2] = grav*cos(u[0]);
		
		dx[0] *= p[0];
		dx[1] *= p[0];
		dx[2] *= p[0];
	}

	bool ode_structure(tw::DiffStructure &s) override {

		s(0,x_indexode(2));
		s(0,u_indexode(0));
		s(0,p_indexode(0));

		s(1,x_indexode(2));
		s(1,u_indexode(0));
		s(1,p_indexode(0));

		s(2,u_indexode(0));
		s(2,p_indexode(0));

		return true;
	}

	
	bool ode_diff(tw::DiffStructure &s, double t, const double *x, const double *u, const double *p) override {

		s(0,x_indexode(2)) = sin(u[0])*p[0];
		s(0,u_indexode(0)) = x[0]*cos(u[0])*p[0];
		s(0,p_indexode(0)) = x[2]*sin(u[0]);

		s(1,x_indexode(2)) = -cos(u[0]);
		s(1,u_indexode(0)) = x[2]*sin(u[0]);
		s(1,p_indexode(0)) = -x[2]*cos(u[0]);

		s(2,u_indexode(0)) = -grav*sin(u[0]);
		s(2,p_indexode(0)) = grav*cos(u[0]);
		
		return true;
	}

	void p_boundary(double *p_low, double *p_upp) override {
		p_low[0] = 0.1;
		p_upp[0] = 5.0;
	}

	void var_boundary(double *x_low, double *x_upp) override {
		
		x_low[x_index(0,0)] = x_upp[x_index(0,0)] = 0;
		x_low[x_index(n_dis-1,0)] = x_upp[x_index(n_dis-1,0)] = 2;
		
		x_low[x_index(0,1)] = x_upp[x_index(0,1)] = 0;
		x_low[x_index(n_dis-1,1)] = x_upp[x_index(n_dis-1,1)] = -2;
		
		x_low[x_index(0,2)] = x_upp[x_index(0,2)] = 0;
	}
};

/////////////////////////////////////////////////////////////////////////////

int main(int argv, char* argc[]) {

	tw::TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv,argc);
	tw::TWfolder folder(&twparameter,0);
	
	folder.worhp_p.TolOpti = 1e-6;
	folder.worhp_p.TolFeas = 2e-10;
	
	tw::TWdimension TWdim;
	TWdim.ID = "Brachistochrone";
	TWdim.n_dis = twparameter.NDIS;
	TWdim.n_ode = 3;
	TWdim.n_ctrl = 1;
	TWdim.n_param = 1;
	TWdim.multinode = {0,TWdim.n_dis/2,TWdim.n_dis-1};
	
	BrachistoPhase ph(TWdim);
	ph.setSolver(&twparameter);
	folder.Add(&ph);

	tw::Viewer *viewer = nullptr;
	if (twparameter.PLOT) viewer = new tw::Viewer(&twparameter);

	folder.Init();
	folder.Init(viewer);
	
	folder.Loop();
	
	delete viewer;

	return 0;
}
