/*----------------------------------------------------------------
 *
 *  Example: Einparken_LKW
 *
 *----------------------------------------------------------------*/

#ifdef WIN32
#include "windows.h"
#endif

#include "TransWORHP.h"

using namespace std;

void einparken2d(glObject *obj, double *x, double t);

double phasenPlot (int &len, int &ndgl, int &nsteuer, double *t,double *x, double *u, int &i, int &index);

class EinparkenPhase : public TransWorhp {
public:

	double b;

	XMLNode* scenenode;

	EinparkenPhase(XMLNode *xmlmain, int dis) : TransWorhp("Einparken mit LKW", dis,7,2,1,0,0) {
		freetime = 1;

		if (xmlmain) {
			scenenode = xmlmain->GetFirstChild("SCENE");
		}
	}


	void OpenWindows(Viewer *viewer) override {
		viewer->PhasePlot("Phasen-Plot", phasenPlot,0,1);
		viewer->ThreeD("Einparken", scenenode, einparken2d);
	}

	string GetXTitle(int d) override {
		if (d==0) return "x";
		if (d==1) return "y";
		if (d==2) return "v";
		if (d==3) return "theta";
		if (d==4) return "phi";
		if (d==5) return "hilfe";
		if (d==6) return "Winkel";
	}

	string GetUTitle(int d) override {
		if (d==0) return "Beschleunigung";
		if (d==1) return "sigma";
	}
	
	
	void p_init(double *p) override {
		p[0] = 2;
	}

	
	double obj() override {
		return p(0) + 4.0*x(n_dis-1, 5);
	}

	bool obj_structure(DiffStructure &s) override {
		s(0,x_index(n_dis-1,5));
		s(0,p_index(0));
		return true;
	}
	
	bool obj_diff(DiffStructure &s) override {
		s(0,x_index(n_dis-1,5)) = 4.0;
		s(0,p_index(0)) = 1;
		return true;
	}
	

	void ode(double *dx, double t, const double *x, const double *u, const double *p) override {
		dx[0] = x[2] * cos(x[3]) * p[0];
		dx[1] = x[2] * sin(x[3]) * p[0];
		dx[2] = u[0] * p[0];
		dx[3] = x[2] * sin(x[4]) * p[0] / b;
		dx[4] = u[1] * p[0];
		dx[5] = u[0]*u[0] + u[1]*u[1];
		dx[6] = -x[2]/1 * sin(x[6]) - 1/1 * dx[3] * cos(x[6]) - dx[3];
	}

	bool ode_structure(DiffStructure &s) override {
		s(0,x_indexode(2));
		s(0,x_indexode(3));
		s(0,p_indexode(0));

		s(1,x_indexode(2));
		s(1,x_indexode(3));
		s(1,p_indexode(0));

		s(2,u_indexode(0));
		s(2,p_indexode(0));

		s(3,x_indexode(2));
		s(3,x_indexode(4));
		s(3,p_indexode(0));

		s(4,u_indexode(1));
		s(4,p_indexode(0));

		s(5,u_indexode(0));
		s(5,u_indexode(1));
		return false;
	}

	bool ode_diff(DiffStructure& s, double t, const double* x, const double* u, const double* p) override {
		s(0,x_indexode(2)) = cos(x[3]) * p[0];
		s(0,x_indexode(3)) = -x[2] * sin(x[3]) * p[0];

		s(1,x_indexode(2)) = sin(x[3]) * p[0];
		s(1,x_indexode(3)) = x[2] * cos(x[3]) * p[0];		

		s(2,u_indexode(0)) = p[0];

		s(3,x_indexode(2)) = sin(x[4]) * p[0] / b;
		s(3,x_indexode(4)) = x[2] * cos(x[4]) * p[0] / b;

		s(4,u_indexode(1)) = p[0];

		s(5,u_indexode(0)) = 2*u[0];
		s(5,u_indexode(1)) = 2*u[1];
		return false;
	}
	
	bool ode_diff_p(DiffStructure& s, double t, const double *x, const double *u, const double *p, int index) override {
		s(0,p_indexode(0)) = x[2] * cos(x[3]);
		s(1,p_indexode(0)) = x[2] * sin(x[3]);
		s(2,p_indexode(0)) = u[0];
		s(3,p_indexode(0)) = x[2] * sin(x[4]) / b;
		s(4,p_indexode(0)) = u[1];
		return false;
	}
	
	void x_boundary(double *x_low, double *x_upp) override {
		x_low[0] = -4;
		x_low[1] = -1;

		x_low[2] = -10;
		x_upp[2] = 10;
		x_low[4] = -M_PI/2; // phi_max
		x_upp[4] = M_PI/2;
		x_low[5] = 0;
		x_low[6] = -M_PI/4; // Winkel_max
		x_upp[6] = M_PI/4;
	}
	
	void u_boundary(double *u_low, double *u_upp) override {
		u_low[0] = -10;
		u_upp[0] = 10;
		u_low[1] = -M_PI/2; // sigma_max
		u_upp[1] = M_PI/2;
	}

	void p_boundary(double *p_low, double *p_upp) override {
		p_low[0] = 1;
		p_upp[0] = 10;
	}

	void var_boundary(double *x_low, double *x_upp) override {
		x_low[x_index(0,0)] = x_upp[x_index(0,0)] = 1;
		x_low[x_index(0,1)] = x_upp[x_index(0,1)] = 1;
		x_low[x_index(0,2)] = x_upp[x_index(0,2)] = 0;
		x_low[x_index(0,3)] = x_upp[x_index(0,3)] = 0;
		x_low[x_index(0,4)] = x_upp[x_index(0,4)] = 0;
		x_low[x_index(0,5)] = x_upp[x_index(0,5)] = 0;
		x_low[x_index(0,6)] = x_upp[x_index(0,6)] = M_PI/4;
		x_low[x_index(n_dis-1,0)] = x_upp[x_index(n_dis-1,0)] = -4;
		x_low[x_index(n_dis-1,1)] = x_upp[x_index(n_dis-1,1)] = -1;
		x_low[x_index(n_dis-1,2)] = x_upp[x_index(n_dis-1,2)] = 0;
		x_low[x_index(n_dis-1,3)] = x_upp[x_index(n_dis-1,3)] = 0;
		x_low[x_index(n_dis-1,4)] = x_upp[x_index(n_dis-1,4)] = 0;
		x_low[x_index(n_dis-1,6)] = x_upp[x_index(n_dis-1,6)] = 0;
	}

};

void einparken2d(glObject *obj, double *x, double t) {

	glLineWidth(2);
	glColor3f(.3,.3,.9);


	//glRotatef(180,0,0,1);


	// Wand
	glBegin(GL_TRIANGLE_FAN);
	glVertex2f(-10, -3);
	glVertex2f(-10, -5);
	glVertex2f(10, -5);
	glVertex2f(10, -3);
	glEnd();
	glBegin(GL_TRIANGLE_FAN);
	glVertex2f(-15, 1);
	glVertex2f(-15, -5);
	glVertex2f(-10, -5);
	glVertex2f(-10, 1);
	glEnd();



	// Auto
	float xx=1.5, yy=1;
	glPushMatrix();
	glTranslatef(x[0],x[1],0);	
	glRotatef(x[3]*180/M_PI,0,0,1); // Rotation

	glBegin(GL_LINE_LOOP);
	glVertex2f(-xx, yy);
	glVertex2f(xx, yy);
	glVertex2f(xx, -yy);
	glVertex2f(-xx, -yy);
	glEnd();

	glPopMatrix();

	//Anhaenger
	glPushMatrix();
	glTranslatef(x[0],x[1],0);
	glRotatef(x[3]*180/M_PI,0,0,1); // Rotation
	glTranslatef(-2.5,0,0);
	glRotatef(x[6]*180/M_PI,0,0,1); // Rotation
	glTranslatef(-2,0,0);

	glBegin(GL_LINE_LOOP);
	glVertex2f(-xx, yy);
	glVertex2f(2*xx, yy);
	glVertex2f(2*xx, -yy);
	glVertex2f(-xx, -yy);
	glEnd();
	glPopMatrix();
}

double phasenPlot (int &len, int &ndgl, int &nsteuer, double *t,double *x, double *u, int &i, int &index) {
	if (index==0) return x[i*(ndgl+nsteuer)+1];
	if (index==1) return x[i*(ndgl+nsteuer)];
}

/////////////////////////////////////////////////////////////////////////////

int main(int argv, char* argc[]) {

	TWparameter twparameter("transworhp.xml");
	std::map<std::string,std::string> a = twparameter.Arguments(argv,argc);
	TWfolder folder(&twparameter,0);

	//b = ToDouble(a["s"]);
	

	XMLNode *xml_laufkatze = TWparameter::ReadParams("laufkatze.xml");

	EinparkenPhase ph(xml_laufkatze, twparameter.NDIS);
	ph.b = 1;
	folder.Add(&ph);

	Viewer *viewer = nullptr;
	if (twparameter.PLOT) viewer = new Viewer(&twparameter);
	
	folder.Init();
	folder.Init(viewer);
	folder.Loop();

	delete viewer;
	
	return 0;
}

