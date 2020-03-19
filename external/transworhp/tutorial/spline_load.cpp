/*----------------------------------------------------------------
 *
 *  Tutorial: Splineproblem mit Initialisierung aus Datei
 *
 *----------------------------------------------------------------*/

#include "TWfolder.h"
#include "TWproblem.h"
#include "TransWORHP.h"

class SplineLoad : public tw::TransWorhpProblem {
public:

	SplineLoad(const tw::TWdimension &TWdata) : TransWorhpProblem(TWdata) {}

	double obj() override {

		return x(n_dis-1,2);
	}

	bool obj_structure(tw::DiffStructure &s) override {

		s(0, x_index(n_dis-1,2));
		return true;
	}

	bool obj_diff(tw::DiffStructure &s) override {

		s(0, x_index(n_dis-1,2)) = 1;
		return true;
	}

	void ode(double *dx, double t, const double *x, const double *u,
			 const double *p) override {

		dx[0] = x[1];
		dx[1] = u[0];
		dx[2] = u[0]*u[0];
	}

	bool ode_structure(tw::DiffStructure &s) override {

		s(0, x_indexode(1)); // dx[0] / dx[1]
		s(1, u_indexode(0)); // dx[1] / du[0]
		s(2, u_indexode(0)); // dx[2] / du[0]
		return true;
	}

	bool ode_diff(tw::DiffStructure &s, double t, const double *x,
			      const double *u, const double *p) override {

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
	twparameter.Arguments(argv, argc);

	tw::TWdimension TWdim;
	TWdim.ID = "Spline / Laden, Speichern";
	TWdim.n_ode = 3;
	TWdim.n_ctrl = 1;
	TWdim.n_dis = twparameter.NDIS;

	tw::Viewer *viewer = nullptr;
	if (twparameter.PLOT) viewer = new tw::Viewer(&twparameter);

	// Startschaetzung in Datei erstellen
	{
		SplineLoad ph(TWdim);
		tw::TWfolder folder(&twparameter, 0);
		ph.setSolver(&twparameter);
		folder.Add(&ph);
		folder.Init();
		folder.Init(viewer);

		folder.Loop(1);
		ph.solver->ToMATLAB("spline11.m");
	}

	if (viewer) viewer->closeAll();

	// Startschaetzung aus Datei lesen (und fuer mehr Diskretisierungspunkte interpolieren)
	{
		//Anzahl diskreter Punkte erhoehen
		TWdim.n_dis = 1001;

		SplineLoad ph(TWdim);
		tw::TWfolder folder(&twparameter, 0);
		ph.setSolver(&twparameter);
		folder.Add(&ph);
		folder.Init();
		folder.Init(viewer);

		ph.solver->FromMATLAB("spline11.m");
		folder.Loop(1);
	}

	delete viewer;

	return 0;
}
