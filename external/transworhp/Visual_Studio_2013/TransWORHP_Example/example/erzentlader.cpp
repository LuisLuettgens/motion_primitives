/*----------------------------------------------------------------
 *
 *  Example: Erzentlader
 *
 *----------------------------------------------------------------*/

#ifdef WIN32
#include "windows.h"
#endif

#include "TransWORHP.h"

using namespace std;

class ErzPhase : public TransWorhp {
public:

	ErzPhase(int dis) : TransWorhp("Erzentlader",dis,4,1,1,0,0) {
		freetime = 1;
	}

	string GetXTitle(int d) override {
		if (d==0) return "x_1";
		if (d==1) return "x_2";
		if (d==2) return "x_3";
		if (d==3) return "x_4";
	}

	string GetUTitle(int d) override {
		if (d==0) return "u";
	}


	void p_init(double *p) override {
		p[0] = 6;
	}
 
	double obj() override {
		return p(0);
	}

	bool obj_structure(DiffStructure &s) override {
		s(0, p_index(0));
		return true;
	}

	bool obj_diff(DiffStructure &s) override {
		s(0,p_index(0)) = 1;
		return true;
	}


	void ode(double *dx, double t, const double *x, const double *u, const double *p) override {
		dx[0] = x[1] * p[0];
		dx[1] = u[0] * p[0];
		dx[2] = x[3] * p[0];
		dx[3] =( -x[2] + u[0]) * p[0];
	}

	bool ode_structure(DiffStructure &s) override {
		s(0, x_indexode(1)); // dx[0] / dx[1]
		s(1, u_indexode(0)); // dx[1] / du[0]
		s(2, x_indexode(3)); // dx[2] / du[0]
		s(3, x_indexode(2)); // dx[2] / du[0]
		s(3, u_indexode(0)); // dx[1] / du[0]
		
		s(0, p_indexode(0)); // dx[1] / du[0]
		s(1, p_indexode(0)); // dx[1] / du[0]
		s(2, p_indexode(0)); // dx[1] / du[0]
		s(3, p_indexode(0)); // dx[1] / du[0]
		return true;
	}

	bool ode_diff(DiffStructure &s, double t, const double *x, const double *u, const double *p) override {
		s(0, x_indexode(1)) = p[0]; // dx[0] / dx[1]
		s(1, u_indexode(0)) = p[0]; // dx[1] / du[0]
		s(2, x_indexode(3)) = p[0]; // dx[2] / du[0]
		s(3, x_indexode(2)) = -p[0]; // dx[2] / du[0]
		s(3, u_indexode(0)) = p[0]; // dx[1] / du[0]
		return true;
	}

	bool ode_diff_p(DiffStructure &s, double t, const double *x, const double *u, const double *p, int index) override {
		s(0, p_indexode(0)) = x[1]; // dx[1] / du[0]
		s(1, p_indexode(0)) = u[0]; // dx[1] / du[0]
		s(2, p_indexode(0)) = x[3]; // dx[1] / du[0]
		s(3, p_indexode(0)) = -x[2]+u[0]; // dx[1] / du[0]
		return true;
	}
	

	void p_boundary(double *p_low, double *p_upp) override {
		p_low[0] = 3;
		p_upp[0] = 10;
	}

	void u_boundary(double *u_low, double *u_upp) override {
		u_low[0] = -1;
		u_upp[0] = +1;
	}

	void x_boundary(double *x_low, double *x_upp) override {
		x_upp[1] = 0.2;
	}

	void var_boundary(double *x_low, double *x_upp) override {
		x_low[x_index ( 0,0 ) ] = x_upp[x_index ( 0,0 ) ] = 0;
		x_low[x_index ( 0,1 ) ] = x_upp[x_index ( 0,1 ) ] = 0;
		x_low[x_index ( 0,2 ) ] = x_upp[x_index ( 0,2 ) ] = 0;
		x_low[x_index ( 0,3 ) ] = x_upp[x_index ( 0,3 ) ] = 0;

		x_low[x_index ( n_dis-1,0 ) ] = x_upp[x_index ( n_dis-1,0 ) ] = 1;
		x_low[x_index ( n_dis-1,1 ) ] = x_upp[x_index ( n_dis-1,1 ) ] = 0;
		x_low[x_index ( n_dis-1,2 ) ] = x_upp[x_index ( n_dis-1,2 ) ] = 0;
		x_low[x_index ( n_dis-1,3 ) ] = x_upp[x_index ( n_dis-1,3 ) ] = 0;
	}

};

/////////////////////////////////////////////////////////////////////////////

int main(int argv, char* argc[]) {

	TWparameter twparameter("transworhp.xml");
	map<string,string> args = twparameter.Arguments(argv,argc);
	TWfolder folder(&twparameter,0);


	ErzPhase ph(twparameter.NDIS);
	folder.Add(&ph);


	Viewer *viewer = nullptr;
	if (twparameter.PLOT) viewer = new Viewer(&twparameter);
	
	folder.Init();
	folder.Init(viewer);
	folder.Loop();


	delete viewer;
	
	return 0;
}

