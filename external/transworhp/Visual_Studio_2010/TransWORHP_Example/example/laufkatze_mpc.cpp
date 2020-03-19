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

double start[] = {0,0,0,0,0,5,0,0,0};
double ziel[] =  {8,0,0,0,0,4,0,0};
//double weight[] = {.015, 0.1, 0.015, 0.1, 0.1, 0.015, .1, .1};
double weight[] = {.015, 0.1, 0.015, 200, 500, 0.15, .1, 0.1};

class LaufkatzeMPC : public TransWorhp {
public:
	double grav;
	XMLNode* scenenode;

	LaufkatzeMPC(XMLNode *xmlmain, int dis) : TransWorhp("Laufkatze", dis,9,2,0,0,0) {

		grav = 9.81;

		if (xmlmain)
			scenenode = xmlmain->GetFirstChild("SCENE");

	}

	void OpenWindows(Viewer *viewer) {

		viewer->ThreeD("Hallenansicht", scenenode, laufkatze3d);

	}
	void GetXTitle(int d, char *s) {
		if (d==0) strcpy(s, "Position X");
		if (d==1) strcpy(s, "Geschwindigkeit X");
		if (d==2) strcpy(s, "Verschiebung Last");
		if (d==3) strcpy(s, "Geschwindigkeit Last");
		if (d==4) strcpy(s, "Beschleunigung X");

		if (d==5) strcpy(s, "Position Y");
		if (d==6) strcpy(s, "Geschwindigkeit Y");
		if (d==7) strcpy(s, "Beschleunigung Y");

		if (d==8) strcpy(s, "u^2");
		
	}

	void GetUTitle(int d, char *s) {
		if (d==0) strcpy(s, "Ruck X");
		if (d==1) strcpy(s, "Ruck Y");
		
	}

	void p_init(double *p) {
		//p[0] = 12;
	}

	void x_init(double *x, int i, int dis) {
		x[5] = 5;
	}
	void u_init(double *u, int i, int dis) {
		u[0] = 0;
		u[1] = 0;
	}


	double obj() {
		
		double ret = 0;
		
		for (int i=0;i<8;i++) {
			double tmp = x(n_dis-1,i)-ziel[i];
			ret += tmp*tmp* weight[i];			
		}
			
		return ret + x(n_dis-1,8);
	}


	bool obj_structure(DiffStructure &s) {
	return false;
	s(0,x_index(n_dis-1,8) );
		return true;

	}

	bool obj_diff(DiffStructure &s) {
		return false;
		s(0,x_index(n_dis-1,8)) = 1 ;
		return true;
	}


	void ode(double *dx, double t, const double *x, const double *u, const double *p) {

		dx[0] = x[1];
		dx[1] = x[4];
		dx[2] = x[3];
		dx[3] = (x[4] - (grav-x[7])*x[2]/x[5]);
		dx[4] = u[0];
		dx[5] = x[6];
		dx[6] = x[7];
		dx[7] = u[1];
		dx[8] = (u[0]*u[0] + u[1]*u[1]);

		//double E_pot = grav * sqrt(x[5]*x[5]-x[2]*x[2]);
		//double E_kin = (x[1]+x[3]) * (x[1]+x[3]) /2;
		
		//dx[8] += (E_pot + E_kin) * .5;
	
	}


	// Optional
	bool ode_structure(DiffStructure &s) {

		return false;
		
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

		s(8,u_indexode(0));
		s(8,u_indexode(1));
		s(8,p_indexode(0));
		
		return true;
	}

	bool ode_diff(DiffStructure& s, double t, const double* x, const double* u, const double* p) {
return false;

//return false;
		// dx[0] = x[1] * p[0];
		s(0,x_indexode(1)) = p[0];
		//s.set(0,p_indexode(0)) = x[1];


		// dx[1] = x[4] * p[0];
		s(1,x_indexode(4)) = p[0];
		//s.use(1,p_indexode(0));

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

		// dx[8] = (u[0]*u[0] + u[1]*u[1])  * p[0];
		s(8,u_indexode(0)) = 2*u[0]*p[0];
		s(8,u_indexode(1)) = 2*u[1]*p[0];

		return true;
	}

	bool ode_diff_p(DiffStructure& s, double t, const double *x, const double *u, const double *p, int index) {
return false;

	
		s(0,p_indexode(0)) = x[1];
		s(1,p_indexode(0)) = x[4];
		s(2,p_indexode(0)) = x[3];
		s(3,p_indexode(0)) = x[4] - (grav-x[7])*x[2]/x[5];
		s(4,p_indexode(0)) = u[0];
		s(5,p_indexode(0)) = x[6];
		s(6,p_indexode(0)) = x[7];
		s(7,p_indexode(0)) = u[1];
		s(8,p_indexode(0)) = u[0]*u[0] + u[1]*u[1];
		
		return true;
		
	}


	void u_boundary(double *u_low, double *u_upp) {

		u_low[0] = -1;
		u_upp[0] =  1;
		u_low[1] = -1;
		u_upp[1] =  1;

	}

	void x_boundary(double *x_low, double *x_upp) {

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

		x_low[5] = 3;//0.5;
		x_upp[5] = 6;//15;

		x_low[6] = -3;
		x_upp[6] = 3;

		x_low[7] = -10;
		x_upp[7] = 10;

		x_low[8] = 0;

	}

	void p_boundary(double *p_low, double *p_upp) {
	
	}

	void var_boundary(double *x_low, double *x_upp) {
	
		for (int i=0;i<9;i++) {
			x_low[x_index(0,i)] = start[i];
			x_upp[x_index(0,i)] = start[i];
		}
		
		/*for (int i=0;i<8;i++) {
			x_low[x_index(n_dis-1,i)] = ziel[i];
			x_upp[x_index(n_dis-1,i)] = ziel[i];
		}*/
		
	//	x_low[x_index(n_dis-1,1)] = x_upp[x_index(n_dis-1,1)] = 0;
		x_low[x_index(n_dis-1,2)] = x_upp[x_index(n_dis-1,2)] = 0;
		x_low[x_index(n_dis-1,3)] = x_upp[x_index(n_dis-1,3)] = 0;
		x_low[x_index(n_dis-1,4)] = x_upp[x_index(n_dis-1,4)] = 0;
	//	x_low[x_index(n_dis-1,6)] = x_upp[x_index(n_dis-1,6)] = 0;
		x_low[x_index(n_dis-1,7)] = x_upp[x_index(n_dis-1,7)] = 0;

		//x_low[x_index(0,8)] = 0;
		//x_upp[x_index(0,8)] = 0;

	}

};


/////////////////////////////////////////////////////////////////////////////



int main(int argv, char* argc[]) {

	TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv,argc);
	

	XMLNode *xml_laufkatze = TWparameter::ReadParams("laufkatze.xml");


	Viewer *viewer = 0;
	if (twparameter.PLOT) viewer = new Viewer(&twparameter);

	
	ofstream mpcfile("mpc.dat");

	for (int K=0;K<60;K++) {
		TWfolder folder(&twparameter,0);
		LaufkatzeMPC ph(xml_laufkatze, twparameter.NDIS);
		folder.Add(&ph);
		
		folder.Init();
		folder.Init(viewer);
		
		for (int i=0;i<ph.n_dis;i++) {
			ph.T[i] = i/(ph.n_dis-1.)*4;
		}
		
		folder.Loop(.5,1);

		int index = ph.n_dis/4-1;
		
		for (int i=0;i<9;i++) {
			start[i] = ph.x(index, i);
			mpcfile << start[i] << " ";
		}
		for (int i=0;i<2;i++) {
			mpcfile << ph.u(index,i) << " ";
		}
		mpcfile << endl;
		
		
		viewer->plots.clear();
	}
	

	
	delete viewer;
	delete xml_laufkatze;

	return 0;
}

