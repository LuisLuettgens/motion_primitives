/*----------------------------------------------------------------
 *
 *  Tutorial: Splineproblem mit Integration als Startschaetzung
 *
 *----------------------------------------------------------------*/

#ifdef WIN32
#include "windows.h"
#endif

#include "TransWORHP.h"

using namespace std;

class SplineRuku : public TransWorhp {
public:

	SplineRuku(int dis) : TransWorhp("Spline / Explizite Integration",dis,3,1,0,0,0) {}

	double obj() {

		return x(n_dis-1,2);
	}

	bool obj_structure(DiffStructure &s) {

		s(0, x_index(n_dis-1,2));
		return true;
	}

	bool obj_diff(DiffStructure &s) {

		s(0, x_index(n_dis-1,2)) = 1;
		return true;
	}

	void ode(double *dx, double t, const double *x, const double *u,
			 const double *p) {

		dx[0] = x[1];
		dx[1] = u[0];
		dx[2] = u[0]*u[0];
	}

	bool ode_structure(DiffStructure &s) {

		s(0, x_indexode(1)); // dx[0] / dx[1]
		s(1, u_indexode(0)); // dx[1] / du[0]
		s(2, u_indexode(0)); // dx[2] / du[0]
		return true;
	}

	bool ode_diff(DiffStructure &s, double t, const double *x,
			      const double *u, const double *p) {

		s(0, x_indexode(1))= 1;
		s(1, u_indexode(0))= 1;
		s(2, u_indexode(0))= 2*u[0];
		return true;
	}

	void u_boundary(double *u_low, double *u_upp) {

		u_low[0] = -6;
		u_upp[0] = +6;
	}

	void x_boundary(double *x_low, double *x_upp) {

		x_low[2] = 0;
	}

	void var_boundary(double *x_low, double *x_upp) {

		x_low[x_index(0,0)] = x_upp[x_index(0,0)] = 0;
		x_low[x_index(0,1)] = x_upp[x_index(0,1)] = 1;
		x_low[x_index(0,2)] = x_upp[x_index(0,2)] = 0;

		x_low[x_index(n_dis-1,0)] = x_upp[x_index(n_dis-1,0)] = 0;
		x_low[x_index(n_dis-1,1)] = x_upp[x_index(n_dis-1,1)] = 1;
	}

};

///////////////////////////////////////////////////////////////////

int main(int argv, char* argc[]) {

	TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv,argc);
	TWfolder folder(&twparameter,0);

	SplineRuku ph(twparameter.NDIS);
	folder.Add(&ph);

	Viewer *viewer = 0;
	if (twparameter.PLOT) viewer = new Viewer(&twparameter);

	folder.Init();
	folder.Init(viewer);

	// Startschaetzung aus Anfangszustand und Steuerung integrieren
	// Anfangszustand festlegen
	ph.X[ph.x_index(0,0)]  = 0;
	ph.X[ph.x_index(0,1)]  = 1;
	ph.X[ph.x_index(0,2)]  = 0;

	// Steuerung festlegen
	for (int i=0;i<twparameter.NDIS;i++) ph.X[ph.u_index(i,0)] = -6 + 12 *  i/(twparameter.NDIS-1.);
		
	// System integrieren
	int steps = ph.Integrate(twparameter.butchertableau);
	cout << "Integrationsschritte: " << steps << endl;

	folder.Loop(5);

	delete viewer;

	return 0;
}
