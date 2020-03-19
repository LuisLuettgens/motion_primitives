/*----------------------------------------------------------------
 *
 * Laufkatze: Einfache Trajektorie
 *
 *----------------------------------------------------------------*/

#ifdef WIN32
#include "windows.h"
#endif

#include "TransWORHP.h"
#include "laufkatze_base.h"

using namespace std;




class Laufkatze_IntPhase : public TransWorhp {
public:
	double grav;
	XMLNode* scenenode;

	Laufkatze_IntPhase(XMLNode *xmlmain, int dis) : TransWorhp("Laufkatze", dis,8,2,1,0,0,1) {
		freetime = 1;
		grav = 9.81;

		if (xmlmain)
			scenenode = xmlmain->GetFirstChild("SCENE");
	}

	void OpenWindows(Viewer *viewer) override {
		viewer->ThreeD("Hallenansicht", scenenode, laufkatze3d);

	}

	string GetXTitle(int d) override {
		if (d==0) return "Position X";
		if (d==1) return "Geschwindigkeit X";
		if (d==2) return "Verschiebung Last";
		if (d==3) return "Geschwindigkeit Last";
		if (d==4) return "Beschleunigung X";

		if (d==5) return "Position Y";
		if (d==6) return "Geschwindigkeit Y";
		if (d==7) return "Beschleunigung Y";

		if (d==8) return "u^2";
	}

	string GetUTitle(int d) override {
		if (d==0) return "Ruck X";
		if (d==1) return "Ruck Y";
	}


	void p_init(double *p) override {
		p[0] = 12;
	}

	void x_init(double *x, int i, int dis) override {
		x[5] = 5;
	}
	void u_init(double *u, int i, int dis) override {
		u[0] = 0;
		u[1] = 0;
	}


	double obj() override {
		return .1 * p(0);
	}

	bool obj_structure(DiffStructure &s) override {
		s(0,p_index(0) );
		return true;
	}

	bool obj_diff(DiffStructure &s) override {
		s(0,p_index(0)) = .1 ;
		return true;
	}


	void ode(double *dx, double t, const double *x, const double *u, const double *p) override {
		dx[0] = x[1] * p[0];
		dx[1] = x[4] * p[0];
		dx[2] = x[3] * p[0];
		dx[3] = (x[4] - (grav-x[7])*x[2]/x[5])  * p[0];
		dx[4] = u[0] * p[0];
		dx[5] = x[6] * p[0];
		dx[6] = x[7] * p[0];
		dx[7] = u[1] * p[0];
	}

	bool ode_structure(DiffStructure &s) override {
		s(0,x_indexode(1));
		s(0,p_indexode(0));

		s(1,x_indexode(4));
		s(1,p_indexode(0));

		s(2,x_indexode(3));
		s(2,p_indexode(0));

		s(3,x_indexode(4));
		s(3,x_indexode(7));
		s(3,x_indexode(2));
		s(3,x_indexode(5));
		s(3,p_indexode(0));

		s(4,u_indexode(0));
		s(4,p_indexode(0));

		s(5,x_indexode(6));
		s(5,p_indexode(0));

		s(6,x_indexode(7));
		s(6,p_indexode(0));

		s(7,u_indexode(1));
		s(7,p_indexode(0));

		return true;

	}

	bool ode_diff(DiffStructure& s, double t, const double* x, const double* u, const double* p) override {
		// dx[0] = x[1] * p[0];
		s(0,x_indexode(1)) = p[0];

		// dx[1] = x[4] * p[0];
		s(1,x_indexode(4)) = p[0];

		// dx[2] = x[3] * p[0];
		s(2,x_indexode(3)) = p[0];

		// dx[3] = (x[4] - (g-x[7])*x[2]/x[5])  * p[0];
		s(3,x_indexode(4)) = p[0];
		s(3,x_indexode(7)) = x[2]/x[5]* p[0];
		s(3,x_indexode(2)) = -(grav-x[7])/x[5]* p[0];
		s(3,x_indexode(5)) = (grav-x[7])*x[2]/x[5]/x[5]* p[0];

		// dx[4] = u[0] * p[0];
		s(4,u_indexode(0)) = p[0];

		// dx[5] = x[6] * p[0];
		s(5,x_indexode(6)) = p[0];

		// dx[6] = x[7] * p[0];
		s(6,x_indexode(7)) = p[0];

		// dx[7] = u[1] * p[0];
		s(7,u_indexode(1)) = p[0];

		return true;
	}

	bool ode_diff_p(DiffStructure& s, double t, const double *x, const double *u, const double *p, int index) override {
		s(0,p_indexode(0)) = x[1];
		s(1,p_indexode(0)) = x[4];
		s(2,p_indexode(0)) = x[3];
		s(3,p_indexode(0)) = x[4] - (grav-x[7])*x[2]/x[5];
		s(4,p_indexode(0)) = u[0];
		s(5,p_indexode(0)) = x[6];
		s(6,p_indexode(0)) = x[7];
		s(7,p_indexode(0)) = u[1];
		return true;	
	}

	void integral(double *f, double t, const double *x, const double *u,
			 const double *p) override {
			
		f[0] = (u[0]*u[0] + u[1]*u[1]) * p[0];
	}

	bool integral_structure(DiffStructure &s) override {
		s(0, u_indexode(0));
		s(0, u_indexode(1));
		s(0, p_indexode(0));
		return true;
	}
	
	bool integral_diff(DiffStructure &s, double t, const double *x, const double *u,
			const double *p) override {
		s(0, u_indexode(0)) = 2*u[0]*p[0];
		s(0, u_indexode(1)) = 2*u[1]*p[0];
		s(0, p_indexode(0)) = u[0]*u[0] + u[1]*u[1];
		return true;
	}
	
	void u_boundary(double *u_low, double *u_upp) override {
		u_low[0] = -1;
		u_upp[0] =  1;
		u_low[1] = -1;
		u_upp[1] =  1;
	}

	void x_boundary(double *x_low, double *x_upp) override {
		x_low[0] = 0;
		x_upp[0] = 100;

		x_low[1] = -3;
		x_upp[1] = 3;

		x_low[2] = -2;
		x_upp[2] = 2;

		x_low[3] = -10;
		x_upp[3] = 10;

		x_low[4] = -4;
		x_upp[4] = 4;

		x_low[5] = 0.5;
		x_upp[5] = 15;

		x_low[6] = -3;
		x_upp[6] = 3;

		x_low[7] = -10;
		x_upp[7] = 10;

	}

	void p_boundary(double *p_low, double *p_upp) override {
		p_low[0] = 1;
		p_upp[0] = 15;
	}

	void var_boundary(double *x_low, double *x_upp) override {

		double start[] = {0,0,0,0,0,5,0,0};
		double ziel[] =  {8,0,0,0,0,4,0,0};

		for (int i=0;i<8;i++) {
			x_low[x_index(0,i)] = start[i];
			x_upp[x_index(0,i)] = start[i];
		}
		for (int i=0;i<8;i++) {
			x_low[x_index(n_dis-1,i)] = ziel[i];
			x_upp[x_index(n_dis-1,i)] = ziel[i];
		}

	}

};


/////////////////////////////////////////////////////////////////////////////


int main(int argv, char* argc[]) {

	TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv,argc);
	TWfolder folder(&twparameter,0);

	XMLNode *xml_laufkatze = TWparameter::ReadParams("laufkatze.xml");

	Laufkatze_IntPhase ph(xml_laufkatze, twparameter.NDIS);
	folder.Add(&ph);

	Viewer *viewer = nullptr;
	if (twparameter.PLOT) viewer = new Viewer(&twparameter);

	folder.Init();
	folder.Init(viewer);

	// Optional: ziemlich gute Startschätzung
	for (int i=0;i<twparameter.NDIS;i++) {
		double z = i/(twparameter.NDIS-1.);
		ph.X[ph.u_index(i,0)] = 6*(z-.5)*(z-.5) - .5;
		ph.X[ph.u_index(i,1)] = -.6*(z-.5)*(z-.5) + .05;
	}
	ph.X[ph.p_index(0)] = 7.9;
	
	ph.Integrate(twparameter.butchertableau);
	// Ende Optional
	
	folder.Loop(2);

	
	delete viewer;
	delete xml_laufkatze;

	return 0;
}

