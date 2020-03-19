/*----------------------------------------------------------------
 *
 *  Example: Raketenwagen mit Gitter-Anpassung
 *
 *----------------------------------------------------------------*/

#ifdef WIN32
#include "windows.h"
#endif

#include "TransWORHP.h"

#include "../src/core/TWparameter.h"
#include "../src/core/Viewer.h"


using namespace std;


class RaketenPhase : public TransWorhp {
public:


	RaketenPhase(int dis) : TransWorhp("Raketenwagen", dis,2,1,1,0,0) {
		freetime = 1;
	}

	void p_init(double *p) override {

		p[0] = 3;

	}

	double obj() override {
		return p(0);
	}


	bool obj_structure(DiffStructure &s) override {

		s(0,p_index(0) );
		return true;

	}

	void ode(double *dx, double t, const double *x, const double *u, const double *p) override {
		dx[0] = x[1] * p[0];
		dx[1] = u[0] * p[0];
	}


	// Optional
	bool ode_structure(DiffStructure &s) override {

		s(0,x_indexode(1));
		s(0,p_indexode(0));

		s(1,u_indexode(0));
		s(1,p_indexode(0));
		
		return true;
	}


	bool ode_diff(DiffStructure& s, double t, const double* x, const double* u, const double* p) override {

		s(0,x_indexode(1)) = p[0];
		s(1,u_indexode(0)) = p[0];
		return true;
	}
	
	bool ode_diff_p(DiffStructure& s, double t, const double *x, const double *u, const double *p, int index) override {
	
		s(0,p_indexode(0)) = x[1];
		s(1,p_indexode(0)) = u[0];
		return true;
	}
		
	void u_boundary(double *u_low, double *u_upp) override {

		u_low[0] = -1;
		u_upp[0] = +1;

	}

	void x_boundary(double *x_low, double *x_upp) override {

	}

	void p_boundary(double *p_low, double *p_upp) override {
		p_low[0] = 1;
		p_upp[0] = 10;
	}

	void var_boundary(double *x_low, double *x_upp) override {

		x_low[x_index(0,0)] = x_upp[x_index(0,0)] = 4;
		x_low[x_index(0,1)] = x_upp[x_index(0,1)] = -1;
		x_low[x_index(n_dis-1,0)] = x_upp[x_index(n_dis-1,0)] = 0;
		x_low[x_index(n_dis-1,1)] = x_upp[x_index(n_dis-1,1)] = 0;
		
	}

};



/////////////////////////////////////////////////////////////////////////////


int main(int argv, char* argc[]) {

	TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv, argc);

	Viewer *viewer = nullptr;
	
	if (twparameter.PLOT) viewer = new Viewer(&twparameter);

	RaketenPhase ph(twparameter.NDIS);

	TWfolder folder(&twparameter, 0);
	folder.Add(&ph);
	folder.Init();
	folder.Init(viewer);
	
	folder.meshRef();

	delete viewer;

	return 0;
}
