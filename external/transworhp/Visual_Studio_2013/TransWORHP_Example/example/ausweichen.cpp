/*----------------------------------------------------------------
 *
 *  Example: Ausweichen
 *
 *----------------------------------------------------------------*/

#ifdef WIN32
#include "windows.h"
#endif

#include "TransWORHP.h"

using namespace std;

void ausweichen2d(glObject *obj, double *x, double t);
double phasenPlot (int &len, int &ndgl, int &nsteuer, double *t,double *x, double *u, int &i, int &index);

class AusweichenPhase : public TransWorhp {
public:

	double Flf;
	double Fsf;
	double Flr;
	double Fsr;

	double FAx;
	double FAy;

	double Jz;
	
	double lr;
	double lf;

	double m;

	double eSP;

	double obsX, obsY;

	double gewichtung1, gewichtung2, gewichtung3;

	double scaleU1;

	XMLNode* scenenode;

	AusweichenPhase(XMLNode *xmlmain, int dis) : TransWorhp("Ausweichen", dis,8,2,1,0,1) {
		freetime = 1;

		if (xmlmain) {
			scenenode = xmlmain->GetFirstChild("SCENE");
		}

		Flf = 1; // Kraefte auf Reifen	(long., vorn)
		Fsf = 1; //         "		(seit., vorn)
		Flr = 1; //         "		(long., hinten)
		Fsr = 1; //         "		(seit., hinten)

		FAx = 1; // Luftwiderstand
		FAy = 1; // Luftwiderstand

		Jz = 1752; //Traegheitsmoment
	
		lr = 1.37484; // Abstand vom COG zum Rad (hinen)
		lf = 1.19016; // Abstand vom COG zum Rad (vorn)

		m = 1300; // Masse

		eSP = 0.5; // Abstand
		
		scaleU1 = 10000;

		// Position Hindernis
		obsX = 20;
		obsY = 1;

	}


	void OpenWindows(Viewer *viewer) override {
		viewer->PhasePlot("Phasen-Plot", phasenPlot,0,1);
		viewer->ThreeD("Ausweichen", scenenode, ausweichen2d);
	}

	string GetXTitle(int d) override {
		if (d==0) return "x_E";
		if (d==1) return "y_E";
		if (d==2) return "v";
		if (d==3) return "phi";
		if (d==4) return "omega_phi";
		if (d==5) return "beta";
		if (d==6) return "delta";
		if (d==7) return "hilfe";
	}

	string GetUTitle(int d) override {
		if (d==0) return "omega_delta";
		if (d==1) return "F_B";
	}
	
	void x_init(double *x, int i, int dis) override {
		x[7] = 0;
	}

	void u_init(double *u, int i, int dis) override {
		//u[0] = 0.6;
		//u[1] = 10000;
	}
	
	void p_init(double *p) override {
		p[0] = 2;
	}

	
	double obj() override {

		gewichtung1 = 2.0;
		gewichtung2 = 5.0;
		gewichtung3 = 0.0;

		return gewichtung1*x(n_dis-1, 0) + gewichtung2*p(0) + gewichtung3*x(n_dis-1, 7);
	}

	bool obj_structure(DiffStructure &s) override {
		s(0,x_index(n_dis-1,0));
		s(0,x_index(n_dis-1,7));
		s(0,p_index(0));
		return true;
	}
	
	bool obj_diff(DiffStructure &s) override {
		s(0,x_index(n_dis-1,0)) = gewichtung1;
		s(0,x_index(n_dis-1,7)) = gewichtung3;
		s(0,p_index(0)) = gewichtung2;
		return true;
	}
	

	void ode(double *dx, double t, const double *x, const double *u, const double *p) override {
		
		double cat1 = 3000;
		double cat2 = 30;
		
		Fsf = cat1*atan(cat2*x[2]);
		Fsr = cat1*atan(cat2*x[2]);

		Flf = (-2.0/3.0*scaleU1*u[1]-m*lf*9.81/(lf+lr));
		Flr = (-1.0/3.0*scaleU1*u[1]-m*lr*9.81/(lf+lr));

		dx[0] = x[2] * cos(x[5]+x[3]) * p[0];
		dx[1] = x[2] * sin(x[5]+x[3]) * p[0];
		dx[2] = ( cos(x[5]) * (Flf*cos(x[6]) - Fsf * sin(x[6]) - FAx + Flr) + sin(x[5]) * (Fsf * cos(x[6]) + Flf*sin(x[6]) + FAy + Fsr) ) * p[0]/m;
		dx[3] = x[4] * p[0];
		dx[4] = ( (Fsf * cos(x[6]) + Flf*sin(x[6]))*lf - Fsf*lr + FAy*eSP ) * p[0]/Jz;
		dx[5] = ( cos(x[5])*(Fsf*cos(x[6]) + Flf*sin(x[6]) + FAy +Fsr) + sin(x[5])*(Fsf*sin(x[6])-Flf*cos(x[6])+FAx-Flr) ) * p[0]/(m*x[2]) - x[4]*p[0];
		dx[6] = u[0] * p[0];
		dx[7] = 0;//(u[0]*u[0] + scaleU1*scaleU1*u[1]*u[1])*p[0];
	}

	bool ode_structure(DiffStructure &s) override {
		/*
		s(0,x_indexode(2));
		s(0,u_indexode(0));
		s(0,p_indexode(0));

		s(1,x_indexode(2));
		s(1,u_indexode(0));
		s(1,p_indexode(0));

		s(2,x_indexode(3));
		s(2,u_indexode(0));
		s(2,p_indexode(0));

		s(3,u_indexode(1));
		s(3,p_indexode(0));

		s(4,u_indexode(0));
		s(4,u_indexode(1));
		*/
		return false;
	}

	bool ode_diff(DiffStructure& s, double t, const double* x, const double* u, const double* p) override {
		/*
		s(0,x_indexode(2)) = -u[0] * sin(x[2]) * p[0];
		s(0,u_indexode(0)) = cos(x[2]) * p[0];

		s(1,x_indexode(2)) = u[0] * cos(x[2]) * p[0];
		s(1,u_indexode(0)) = sin(x[2]) * p[0];

		s(2,x_indexode(3)) = u[0] * cos(x[3]) * p[0];
		s(2,u_indexode(0)) = sin(x[3]) * p[0] / 1.0;

		s(3,u_indexode(1)) = p[0];

		s(4,u_indexode(0)) = 2.0*u[0];
		s(4,u_indexode(1)) = 2.0*u[1];
		*/
		return false;
	}
	
	bool ode_diff_p(DiffStructure& s, double t, const double *x, const double *u, const double *p, int index) override {
		/*		
		s(0,p_indexode(0)) = u[0] * cos(x[2]);
		s(1,p_indexode(0)) = u[0] * sin(x[2]);
		s(2,p_indexode(0)) = u[0] * sin(x[3]) / 1.0;
		s(3,p_indexode(0)) = u[1];
		*/
		return false;
	}
	
	void x_boundary(double *x_low, double *x_upp) override {
		
		x_low[0] = 0;

		x_low[1] = -10;
		x_upp[1] = 10;

		x_low[2] = 10;
		x_upp[2] = 100;

		x_low[3] = -M_PI/2.0;
		x_upp[3] =  M_PI/2.0;

		x_low[5] = -M_PI/2.0;
		x_upp[5] =  M_PI/2.0;

		x_low[6] = -M_PI/2.0;
		x_upp[6] =  M_PI/2.0;

		x_low[7] = 0;
	}
	
	void u_boundary(double *u_low, double *u_upp) override {
		u_low[0] = -M_PI/4.;
		u_upp[0] =  M_PI/4.;

		u_low[1] = 0.0;
		u_upp[1] = 15000/scaleU1;
	}

	void p_boundary(double *p_low, double *p_upp) override {
		p_low[0] = 0.10;
		p_upp[0] = 10.0;
	}

	void var_boundary(double *x_low, double *x_upp) override {
		x_low[x_index(0,0)] = x_upp[x_index(0,0)] = 0;
		x_low[x_index(0,1)] = x_upp[x_index(0,1)] = 0;
		x_low[x_index(0,2)] = x_upp[x_index(0,2)] = 70;
		x_low[x_index(0,3)] = x_upp[x_index(0,3)] = 0;
		x_low[x_index(0,4)] = x_upp[x_index(0,4)] = 0;
		x_low[x_index(0,5)] = x_upp[x_index(0,5)] = 0;
		x_low[x_index(0,6)] = x_upp[x_index(0,6)] = 0;
		
		x_low[x_index(n_dis-1,0)] = x_upp[x_index(n_dis-1,0)] = 20;
		//x_low[x_index(n_dis-1,1)] = x_upp[x_index(n_dis-1,1)] = 0;
		x_low[x_index(n_dis-1,3)] = x_upp[x_index(n_dis-1,3)] = 0;
		x_low[x_index(n_dis-1,5)] = x_upp[x_index(n_dis-1,5)] = 0;
		x_low[x_index(n_dis-1,6)] = x_upp[x_index(n_dis-1,6)] = 0;
	}
	
	void rand(double *r) override {
		r[0] = x(n_dis-1,3) + x(n_dis-1,5);
		r[1] = obsX - x(n_dis-1,0);
	}
	
	void rand_boundary(double *r_low, double *r_upp) override {
		r_upp[0] = 0;
		r_upp[1] = 0;
	}
	
	void neben(double *c, double t, const double *x, const double *u, const double *p) override {
		c[0] = ((x[0]-obsX)*(x[0]-obsX) + (x[1]-obsY)*(x[1]-obsY));
	}
	
	void neben_boundary(double *c_low, double *c_upp) override {
		double rE = 1.5;
		double rObs = 1.5;
		c_low[0] = (rE + rObs)*(rE + rObs);
	}
};

void ausweichen2d(glObject *obj, double *x, double t) {

	glLineWidth(2);
	glColor3f(.3,.3,.9);

	// Auto
	float xx=1.5, yy=1;
	
	glTranslatef(x[0],x[1],0);	
	glRotatef(x[5]*180/M_PI,0,0,1); // Rotation

	glBegin(GL_LINE_LOOP);
	glVertex2f(-xx, yy);
	glVertex2f(xx, yy);
	glVertex2f(xx, -yy);
	glVertex2f(-xx, -yy);
	glEnd();

	
}

double phasenPlot (int &len, int &ndgl, int &nsteuer, double *t,double *x, double *u, int &i, int &index) {
	if (index==0) return x[i*(ndgl+nsteuer)+1];
	if (index==1) return x[i*(ndgl+nsteuer)];
}

/////////////////////////////////////////////////////////////////////////////

int main(int argv, char* argc[]) {

	TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv,argc);
	TWfolder folder(&twparameter,0);

	XMLNode *xml_laufkatze = TWparameter::ReadParams("laufkatze.xml");

	AusweichenPhase ph(xml_laufkatze, twparameter.NDIS);
	folder.Add(&ph);

	Viewer *viewer = nullptr;
	if (twparameter.PLOT) viewer = new Viewer(&twparameter);

	folder.Init();
	folder.Init(viewer);

	
	ph.X[ph.x_index(0,2)] = 70;
	for (int i=0;i<twparameter.NDIS;i++) {
		double z = i/(twparameter.NDIS-1.);
		ph.X[ph.u_index(i,0)] = (i>twparameter.NDIS/4 && i<twparameter.NDIS*0.75/4) ? -0.7 : 0.7 ;//6*(z-.5)*(z-.5) - .5;
		ph.X[ph.u_index(i,1)] = 0;//-.6*(z-.5)*(z-.5) + .05;
	}
	ph.X[ph.p_index(0)] = 1;
	
	ph.Integrate(twparameter.butchertableau);
	

	folder.Loop();

	delete viewer;
	
	return 0;
}

