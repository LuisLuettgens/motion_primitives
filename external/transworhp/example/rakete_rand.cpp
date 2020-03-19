/*----------------------------------------------------------------
 *
 *  Example: Raketenwagen mit Randbedingungen
 *
 *----------------------------------------------------------------*/

#include "TWfolder.h"
#include "TWproblem.h"
#include "TransWORHP.h"

class RaketenPhase : public tw::TransWorhpProblem {
public:

	RaketenPhase(const tw::TWdimension &TWdata) : TransWorhpProblem(TWdata) {
		freetime = true;
	}

	void selectWindows(tw::Viewer *viewer) override {

		viewer->AddStateView(0,"x_1");
		viewer->AddStateView(1,"x_2");

		viewer->AddControlView(0,"u");
	}

	void p_init(double *p) override {
		p[0] = 3.0;
	}

	void u_init(double *u, int i, int dis) override {

		if (i < dis*0.5) {
			u[0] = -1.0;
		} else {
			u[0] = 1.0;
		}
	}

	double obj() override {
		return p(0);
	}

	bool obj_structure(tw::DiffStructure &s) override {
		s(0,p_index(0));
		return true;
	}

	bool obj_diff(tw::DiffStructure &s) override {
		s(0,p_index(0)) = 1.0;
		return true;
	}


	void ode(double *dx, double t, const double *x, const double *u, const double *p) override {
		dx[0] = x[1] * p[0];
		dx[1] = u[0] * p[0];
	}

	bool ode_structure(tw::DiffStructure &s) override {
		s(0,x_indexode(1));
		s(0,p_indexode(0));

		s(1,u_indexode(0));
		s(1,p_indexode(0));

		return true;
	}

	bool ode_diff(tw::DiffStructure& s, double t, const double* x, const double* u, const double* p) override {
		s(0,x_indexode(1)) = p[0];
		s(1,u_indexode(0)) = p[0];

		return true;
	}

	bool ode_diff_p(tw::DiffStructure& s, double t, const double *x, const double *u, const double *p, int index) override {
		s(0,p_indexode(0)) = x[1];
		s(1,p_indexode(0)) = u[0];

		return true;
	}


	void u_boundary(double *u_low, double *u_upp) override {
		u_low[0] = -1.0;
		u_upp[0] = +1.0;
	}


	void rand(double *r) override {
		r[0] = x(0,0);
		r[1] = x(0,1);
		r[2] = x(n_dis-1,0);
		r[3] = x(n_dis-1,1);
		r[4] = p(0);
	}

	void rand_boundary(double *r_low, double *r_upp) override {
		r_low[0] = r_upp[0] = 4.0;
		r_low[1] = r_upp[1] = -1.0;
		r_low[2] = r_upp[2] = 0.0;
		r_low[3] = r_upp[3] = 0.0;
		r_low[4] = 1.0;
		r_upp[4] = 10.0;
	}

	bool rand_structure(tw::DiffStructure &s) override {
		s(0,x_index(0,0));
		s(1,x_index(0,1));
		s(2,x_index(n_dis-1,0));
		s(3,x_index(n_dis-1,1));
		s(4,p_index(0));

		return true;
	}

	bool rand_diff(tw::DiffStructure &s) override {
		s(0,x_index(0,0)) = 1.0;
		s(1,x_index(0,1)) = 1.0;
		s(2,x_index(n_dis-1,0)) = 1.0;
		s(3,x_index(n_dis-1,1)) = 1.0;
		s(4,p_index(0)) = 1.0;

		return true;
	}
};

/////////////////////////////////////////////////////////////////////////////

int main(int argv, char* argc[]) {

	tw::TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv,argc);
	tw::TWfolder folder(&twparameter,0);

	tw::TWdimension TWdim;
	TWdim.ID = "Raketenwagen mit Randbedingungen";
	TWdim.n_dis = twparameter.NDIS;
	TWdim.n_ode = 2;
	TWdim.n_ctrl = 1;
	TWdim.n_param = 1;
	TWdim.n_rand = 5;
	TWdim.setMultinodes(3);

	RaketenPhase ph(TWdim);
	ph.setSolver(&twparameter);
	folder.Add(&ph);

	std::unique_ptr<tw::Viewer> viewer;
	if (twparameter.PLOT) viewer = std::unique_ptr<tw::Viewer>(new tw::Viewer(&twparameter));

	folder.Init();
	folder.Init(viewer.get());

	folder.Loop();

	return 0;
}
