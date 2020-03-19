/*----------------------------------------------------------------
 *
 * Laufkatze: Einfache Trajektorie
 *
 *----------------------------------------------------------------*/

#ifdef WIN32
#include "windows.h"
#endif

#include "TransWORHP.h"

using namespace std;


void laufkatze3d(glObject *obj, double *x, double t) {

	float offset = 4;

	int n = obj->countObj();


	Vektor<float> v(0,0,0);

	obj->Draw(v,0,0,1); // Ebene
	obj->Draw(v,0,0,2); // Regal1
	obj->Draw(v,0,0,3); // Regal2


	glPushMatrix();
	glTranslatef(x[0], 0, 0);
	obj->Draw(v,0,0,4);
	obj->Draw(v,0,0,5);

	obj->Draw(v,0,0,8);
	obj->Draw(v,0,0,9);
	obj->Draw(v,0,0,10);
	obj->Draw(v,0,0,11);

	glPopMatrix();

	glPushMatrix();
	glTranslatef(x[0]-x[2], 0, 9 - x[5]);
	obj->Draw(v,0,0,6);
	obj->Draw(v,0,0,7);
	glPopMatrix();


	// Seile
	glLineWidth(2);
	glColor3f(.3,.3,.9);
	glDisable(GL_LIGHTING);

	float xx=1.2, yy=1;

	glBegin(GL_LINES);
	glVertex3f(x[0] +xx, +yy, 9.2);
	glVertex3f(x[0]-x[2]+xx, +yy, 9-x[5]);

	glVertex3f(x[0] -xx, +yy, 9.2);
	glVertex3f(x[0]-x[2]-xx, +yy, 9-x[5]);

	glVertex3f(x[0] +xx, -yy, 9.2);
	glVertex3f(x[0]-x[2]+xx, -yy, 9-x[5]);

	glVertex3f(x[0] -xx, -yy, 9.2);
	glVertex3f(x[0]-x[2]-xx, -yy, 9-x[5]);
	glEnd();

}


class LaufkatzeScalePhase : public ScaleTransWorhp {
public:
	double grav;
	XMLNode* scenenode;
	
	double usc_u[2];
	double usc_l[2];
	double xsc_u[9];
	double xsc_l[9];

	LaufkatzeScalePhase(XMLNode *xmlmain, int dis) : ScaleTransWorhp("Laufkatze Scale", dis,9,2,1,0,0) {
		freetime = 1;
		grav = 9.81;

		
	for (int i=0;i<n_ode;i++) {
		xsc_l[i] = 0;
		xsc_u[i] = 1;
	}
	
	
	usc_u[0] = 	.2;
	usc_u[1] = 	.04;
	usc_l[0] = 	-.3;
	usc_l[1] = -.06;
	
	
	
		if (xmlmain)
			scenenode = xmlmain->GetFirstChild("SCENE");

	}

	void OpenWindows(Viewer *viewer) override {

		viewer->ThreeD("Hallenansicht", scenenode, laufkatze3d);

	}

	void selectWindows(Viewer *viewer) override {

		viewer->AddStateView(0,"Position X");
		viewer->AddStateView(1,"Geschwindigkeit X");
		viewer->AddStateView(2,"Verschiebung Last");
		viewer->AddStateView(3,"Geschwindigkeit Last");
		viewer->AddStateView(4,"Beschleunigung X");
		viewer->AddStateView(5,"Position Y");
		viewer->AddStateView(6,"Geschwindigkeit Y");
		viewer->AddStateView(7,"Beschleunigung Y");
		viewer->AddStateView(8,"u^2");
		
		viewer->AddControlView(0,"Ruck X");
		viewer->AddControlView(1,"Ruck Y");
	}

	void p_init(double *p) override {
		p[0] = 12;
	}

	void SC_x_init(double *x, int i, int dis) override {
		x[5] = 5;
	}
	void SC_u_init(double *u, int i, int dis) override {
		u[0] = 0;
		u[1] = 0;
	}


	double obj() override {
		return x(n_dis-1,8) + .1 * p(0);
	}


	bool obj_structure(DiffStructure &s) override {
		s(0,x_index(n_dis-1,8) );
		s(0,p_index(0) );
		return true;

	}

	bool obj_diff(DiffStructure &s) override {
		s(0,x_index(n_dis-1,8)) = 1 ;
		s(0,p_index(0)) = .1 ;
		return true;
	}


// Check, ob zusammenpasst???


void u_scale(double *u) override {

	for (int i=0;i<n_ctrl;i++)
		u[i] = u[i] * (usc_u[i]-usc_l[i]) + usc_l[i];
	
}

void u_unscale(double *u) override {

	for (int i=0;i<n_ctrl;i++)
		u[i] = (u[i] - usc_l[i]) / (usc_u[i]-usc_l[i]);

}

void x_scale(double *x) override {

	for (int i=0;i<n_ode;i++)
		x[i] = x[i] * (xsc_u[i]-xsc_l[i]) + xsc_l[i];
	
}

void x_unscale(double *x) override {

	for (int i=0;i<n_ode;i++)
		x[i] = (x[i] - xsc_l[i]) / (xsc_u[i]-xsc_l[i]);

}

// bei Plot zurückmachen.



	void SC_ode(double *dx, double t, const double *x, const double *u, const double *p) override {

		dx[0] = x[1] * p[0];
		dx[1] = x[4] * p[0];
		dx[2] = x[3] * p[0];
		dx[3] = (x[4] - (grav-x[7])*x[2]/x[5])  * p[0];
		dx[4] = u[0] * p[0];
		dx[5] = x[6] * p[0];
		dx[6] = x[7] * p[0];
		dx[7] = u[1] * p[0];
		dx[8] = (u[0]*u[0] + u[1]*u[1])  * p[0];

	}


	// Optional
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

		s(8,u_indexode(0));
		s(8,u_indexode(1));
		s(8,p_indexode(0));

		return true;
	}

	bool SC_ode_diff(DiffStructure& s, double t, const double* x, const double* u, const double* p) override {

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

	bool SC_ode_diff_p(DiffStructure& s, double t, const double *x, const double *u, const double *p, int index) override {


	/*	s.set(0,p_indexode(0), x[1]);
		s.set(1,p_indexode(0), x[4]);
		s.set(2,p_indexode(0), x[3]);
		s.set(3,p_indexode(0), x[4] - (grav-x[7])*x[2]/x[5] );
		s.set(4,p_indexode(0), u[0]);
		s.set(5,p_indexode(0), x[x6]);
		s.set(6,p_indexode(0), x[7]);
		s.set(7,p_indexode(0), u[1]);
		s.set(8,p_indexode(0), u[0]*u[0] + u[1]*u[1]);
*/
	
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


	void SC_u_boundary(double *u_low, double *u_upp) override {
		
		u_low[0] = -1;
		u_upp[0] =  1;
		u_low[1] = -1;
		u_upp[1] =  1;

	}

	void SC_x_boundary(double *x_low, double *x_upp) override {

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

		x_low[8] = 0;

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

		x_low[x_index(0,8)] = 0;
		x_upp[x_index(0,8)] = 0;

	}

};


/////////////////////////////////////////////////////////////////////////////



int main(int argv, char* argc[]) {

	TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv,argc);
	TWfolder folder(&twparameter,0);

	XMLNode *xml_laufkatze = TWparameter::ReadParams("laufkatze.xml");



	LaufkatzeScalePhase ph(xml_laufkatze, twparameter.NDIS);
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

