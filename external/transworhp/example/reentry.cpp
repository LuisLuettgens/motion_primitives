/*----------------------------------------------------------------
 *
 * Example: Reentry-Problem
 *
 *----------------------------------------------------------------*/

#include "TWfolder.h"
#include "TWproblem.h"
#include "TransWORHP.h"

#include "twstatus.h"

class ReentryPhase : public tw::TransWorhpProblem {
public:

	double R, beta, p0, g, S_M;

	ReentryPhase(int dis) : TransWorhpProblem(tw::TWdimension("Reentry-Problem",dis,4,1,1,0,0)) {
		freetime = true;

		R = 209;
		beta = 4.26;
		p0 = 2.704e-3;
		g = 3.2172e-4;
		S_M = 53200;
	}

	std::string GetXTitle(int d) override {
		if (d==0) return "v";
		if (d==1) return "gamma";
		if (d==2) return "xi";
		if (d==3) return "integral";
	}

	std::string GetUTitle(int d) override {
		if (d==0) return "u";
	}

	void p_init(double *p) override {
		p[0] = 220;
	}

	void u_init(double *u, int i, int dis) override {
		u[0] = 0.5;
	}

	void x_init(double *x, int i, int dis) override {
		x[0] = 0.36;
		x[1] = -.1;
		x[2] = 0.02;
	}


	double obj() override {
		return x(n_dis-1,3);
	}

	bool obj_structure(tw::DiffStructure &s) override {
		s(0,x_index(n_dis-1,3) );
		return true;
	}

	bool obj_diff(tw::DiffStructure &s) override {
		s(0,x_index(n_dis-1,3) ) = 1;
		return true;
	}


	void ode(double *dx, double t, const double *x, const double *u, const double *p) override {

		double v = x[0];
		double gamma = x[1];
		double xi = x[2];

		double P = p0 * exp(-beta * R * xi);
		double Cw = 1.174-0.9*cos(u[0]);
		double Ca = 0.6*sin(u[0]);


		dx[0] = -S_M * P * v*v / 2 *Cw -g*sin(gamma)/((1+xi)*(1+xi));
		dx[1] = S_M * P * v / 2 *Ca + v*cos(gamma)/(R*(1+xi)) - g*cos(gamma)/(v*(1+xi)*(1+xi));
		dx[2] = v*sin(gamma)/R;
		dx[3] = 10*v*v*v*sqrt(P);

		for (int i=0;i<4;i++)
			dx[i] *= p[0];

	}

	bool ode_structure(tw::DiffStructure &s) override {
		return false;
// 		s(0,u_index(0,0));
// 		s(0,x_index(0,0));
// 		s(0,x_index(0,1));
//
// 		s(1,u_index(0,0));
// 		s(1,x_index(0,0));
// 		s(1,x_index(0,1));
// 		return true;
	}

	bool ode_diff(tw::DiffStructure& s, double t, const double* x, const double* u, const double* p) override {
		return false;
// 		s(0,u_index(0,0)) = -x[0] + 2*u[0]*x[1];
// 		s(0,x_index(0,0)) = -u[0];
// 		s(0,x_index(0,1)) = u[0]*u[0];
//
// 		s(1,u_index(0,0)) = x[0] - 6*u[0]*x[1];
// 		s(1,x_index(0,0)) = u[0];
// 		s(1,x_index(0,1)) = -3*u[0]*u[0];
// 		return true;
	}


	void p_boundary(double *p_low, double *p_upp) override {
		p_low[0] = 200;
	}

	void u_boundary(double *u_low, double *u_upp) override {
		u_low[0] = -10.5;
		u_upp[0] = +10.5;
	}

	void var_boundary(double *x_low, double *x_upp) override {

		double R = 209;

		x_upp[x_index(0,0)] = x_low[x_index(0,0)] = 0.36;
		x_upp[x_index(0,1)] = x_low[x_index(0,1)] = -8.1*M_PI/180;
		x_upp[x_index(0,2)] = x_low[x_index(0,2)] = 4/R;
		x_upp[x_index(0,3)] = x_low[x_index(0,3)] = 0;

		x_upp[x_index(n_dis-1,0)] = x_low[x_index(n_dis-1,0)] = 0.27;
		x_upp[x_index(n_dis-1,1)] = x_low[x_index(n_dis-1,1)] = 0;
		x_upp[x_index(n_dis-1,2)] = x_low[x_index(n_dis-1,2)] = 2.5/R;

	}

	void terminate() override {
		MyStatus("", "Final Time: " + std::to_string(p(0)), Status::NORMAL);
	}


};


class MyFolder : public tw::TWfolder {

	public:

		MyFolder(tw::TWparameter *twp, int a) : TWfolder(twp,a) {}
		bool step() override {

			for (int i=0;i<phases[0]->n_dis;i++) {
				while (phases[0]->u(i,0)>3*M_PI/2) phases[0]->solver->X[phases[0]->u_index(i,0)] -= 2*M_PI;
				while (phases[0]->u(i,0)<-M_PI/2) phases[0]->solver->X[phases[0]->u_index(i,0)] += 2*M_PI;
			}
			return true;
		}
};

/////////////////////////////////////////////////////////////////////////////

int main(int argv, char* argc[]) {

	tw::TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv,argc);
	MyFolder folder(&twparameter,0);

	ReentryPhase ph(twparameter.NDIS);
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
