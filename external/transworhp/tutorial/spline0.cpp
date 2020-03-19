/*----------------------------------------------------------------
 *
 *  Tutorial: Splineproblem ohne Grafik
 *
 *----------------------------------------------------------------*/

#include "TWfolder.h"
#include "TWproblem.h"
#include "TransWORHP.h"

class Spline0 : public tw::TransWorhpProblem {
public:

	Spline0(const tw::TWdimension &TWdata) : TransWorhpProblem(TWdata) {}

	double obj() override {

		return x(n_dis-1,2);
	}

	void ode(double *dx, double t, const double *x, const double *u,
			 const double *p) override {

		dx[0] = x[1];
		dx[1] = u[0];
		dx[2] = u[0]*u[0];
	}

	void u_boundary(double *u_low, double *u_upp) override {

		u_low[0] = -6;
		u_upp[0] = +6;
	}

	void x_boundary(double *x_low, double *x_upp) override {

		x_low[2] = 0;
	}

	void var_boundary(double *x_low, double *x_upp) override {

		x_low[x_index(0,0)] = x_upp[x_index(0,0)] = 0;
		x_low[x_index(0,1)] = x_upp[x_index(0,1)] = 1;
		x_low[x_index(0,2)] = x_upp[x_index(0,2)] = 0;

		x_low[x_index(n_dis-1,0)] = x_upp[x_index(n_dis-1,0)] = 0;
		x_low[x_index(n_dis-1,1)] = x_upp[x_index(n_dis-1,1)] = 1;
	}

};

///////////////////////////////////////////////////////////////////

int main(int argv, char* argc[]) {

	tw::TWparameter twparameter("transworhp.xml");
	tw::TWfolder folder(&twparameter, 0);

	tw::TWdimension TWdim;
	TWdim.ID = "Spline0";
	TWdim.n_ode = 3;
	TWdim.n_ctrl = 1;
	TWdim.n_dis = 11;

	Spline0 ph(TWdim);
	ph.setSolver(&twparameter);
	folder.Add(&ph);

	folder.Init();

	folder.Loop();

	return 0;
}
