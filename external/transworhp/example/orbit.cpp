/*----------------------------------------------------------------
 *
 * GTOC8
 *
 *----------------------------------------------------------------*/

#define _USE_MATH_DEFINES
#include <cmath>

#ifdef WIN32
#include "windows.h"
#endif

#include "TransWORHP.h"
#include "ExplTransWORHP.h"
using namespace std;



double start[] = {4000,0,0,0,8,0,4000};
		
double uscale =1;
/** globaler Zeiger auf OCP für Grafik */
TransWorhp *pph=nullptr;

double flugzeugplot (int &len, int &ndgl, int &nsteuer, double *t,double *x, double *u, int &i, int &index) {

	if (index==0) {
		int g_index = (pph->n_dis-1) * pph->n_ode + pph->n_rand + i*pph->n_neben;
		return pph->G[g_index + 0];
	}

	return 0;
}


void flugzeug3d(glObject *obj, double *x, double t) {

	double scale = 100;
	
	glLineWidth(2);
	glColor3f(.8,.1,.1);
	glDisable(GL_LIGHTING);
	glBegin(GL_LINE_STRIP);
	
	
	for (int i=0;i<pph->n_dis;i++)
		glVertex3f((pph->x(i,0)-400)/scale,pph->x(i,1)/scale,pph->x(i,2)/scale*30 );
	
	glEnd();
	
	
	glEnable(GL_LIGHTING);
	
	float offset = 4;

	int n = obj->countObj();

	Vektor<float> v(0,0,0);
	
	glColor3f(153/255.,185/255.,131/255.);
	obj->Draw(v,0,0,2);
	
	obj->Draw(v,0,0,1);
	
	
	glPushMatrix();
	glTranslatef((x[0]-400)/scale, x[1]/scale, x[2]/scale*30);
	
	glRotatef(x[3]*180/M_PI,0,0,1);
	glRotatef(10* x[4]*180/M_PI,0,1,0);
	glRotatef(-90,0,0,1);
	glColor3f(1,1,1);
	
	obj->Draw(v,0,0,4);
	glPopMatrix();
	/*

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
	glEnd();*/


}


class OrbitPhase : public ExplTransWorhp {
public:
	
	OrbitPhase(int dis, vector<int> multinode) : ExplTransWorhp("Orbit", dis,multinode,7,3,1,0,1) {
		/** freie Endzeit */
		freetime = 1;	

	}

	void OpenWindows(Viewer *viewer) override {

		//viewer->ThreeD("Flugbahn", scenenode, flugzeug3d);
		//viewer->Data("Lastfaktor", flugzeugplot,0);

	}
	/*void GetXTitle(int d, char *s) {
		if (d==0) strcpy(s, "Position X");
		if (d==1) strcpy(s, "Position Y");
		if (d==2) strcpy(s, "Position Z");
		if (d==3) strcpy(s, "Kurswinkel chi");
		if (d==4) strcpy(s, "Steigwinkel gamma");
		if (d==5) strcpy(s, "Haengewinkel mu");
		if (d==6) strcpy(s, "Geschwindigkeit V");
		if (d==7) strcpy(s, "Schubhebelstellung delta_t");
	}

	void GetUTitle(int d, char *s) {
		if (d==0) strcpy(s, "Anstellwinkel alpha_CMD");
		if (d==1) strcpy(s, "Schiebewinkel beta_CMD");
		if (d==2) strcpy(s, "Haengewinkelrate dot_mu_CMD");
		if (d==3) strcpy(s, "Schubhebelvorgabe delta_t_CMD");

	}*/

	void p_init(double *p) override {
		/** Startschaetzung der Flugzeit */ 
		p[0] = 1000;
	}

	void x_init(double *x, int i, int dis) override {
		
		for (int j=0;j<7;j++) {
			x[j] = start[j];
		}
	}
	void u_init(double *u, int i, int dis) override {
		u[0] = 0.5/uscale;
		// u[1] = 0;
	}


	double obj() override {
		/** Mayer-Anteil der Zielfunktion */
		return //1 * p(0) - x(n_dis-1,6) + 
			x(n_dis-1,0)  ;
	}


	bool obj_structure(DiffStructure &s) override {
	
	//s(0,p_index(0));
	//	s(0,x_index(n_dis-1,6));
		s(0,x_index(n_dis-1,0));
		return true;
	}

	bool obj_diff(DiffStructure &s) override {
		s(0,x_index(n_dis-1,0))=1;
		return true;
		s(0,p_index(0)) = 1;
		s(0,x_index(n_dis-1,6)) = -1;
		return true;
	}





	void ode(double *dx, double t, const double *x, const double *u, const double *p) override {

		double mu = 398600.4329;
		double r = sqrt(x[0]*x[0]+ x[1]*x[1]+ x[2]*x[2]);

		double Isp=5000;
		double g=9.80665;
		double m = x[6];

		Vektor<double> VEL(x[3],x[4],x[5]);
		//%VEL = [0,0,1]; % x(4:6)
		Vektor<double> POS(x[0], x[1],x[2]);
		Vektor<double> h = POS.cross(VEL);
		Vektor<double> vv = POS.cross(h);
		double vvnorm = vv.Length();
		double hnorm = h.Length();
		 
		Vektor<double> T;
		T.X() =  uscale *(u[0] * x[0]/r + u[1] * h.X()/hnorm+ u[2] * vv.X()/vvnorm ) ;
		T.Y() =  uscale*(u[0] * x[1]/r + u[1] * h.Y()/hnorm + u[2] * vv.Y()/vvnorm);
		T.Z() =  uscale*(u[0] * x[2]/r + u[1] * h.Z()/hnorm+ u[2] * vv.Z()/vvnorm  );
		  
		dx[0] = x[3];
		dx[1] = x[4];
		dx[2] = x[5];
		dx[3] = -mu * x[0] / (r*r*r) + T.X()/m;
		dx[4] = -mu * x[1] / (r*r*r) + T.Y()/m;
		dx[5] = -mu * x[2] / (r*r*r) + T.Z()/m;
		dx[6] = -uscale * sqrt(u[0]*u[0]+u[1]*u[1]+u[2]+u[2]) / (Isp*g);

		for (int i=0;i<7;i++) {
			dx[i] *= p[0];
		}
	}


	// Optional
	bool ode_structure(DiffStructure &s) override {

		return false;


		for (int i=0;i<8;i++) {

			s(i,p_indexode(0));
		}

		for (int i=1;i<8;i++) {
			if (i==7) continue;
			if (i==5) continue;
			if (i==2) continue;
			if (i==4) continue;
		
			s(i,x_indexode(4));
			s(i,x_indexode(5));
			s(i,x_indexode(6));
			s(i,x_indexode(7));
			s(i,u_indexode(0));
			s(i,u_indexode(1));
			s(i,u_indexode(2));
			s(i,u_indexode(3));
		}

		s(0,x_indexode(3));
		s(0,x_indexode(4));
		s(0,x_indexode(6));
		s(0,p_indexode(0));

		s(1,x_indexode(3));
		s(1,x_indexode(4));
		s(1,x_indexode(6));
		s(1,p_indexode(0));

		s(2,x_indexode(4));
		s(2,x_indexode(6));
		s(2,p_indexode(0));

		s(3,x_indexode(4));
		s(3,x_indexode(5));
		s(3,x_indexode(6));
		s(3,u_indexode(0));
		s(3,u_indexode(1));
		s(3,p_indexode(0));

		s(4,x_indexode(4));
		s(4,x_indexode(5));
		s(4,x_indexode(6));
		s(4,u_indexode(0));
		s(4,u_indexode(1));
		s(4,p_indexode(0));

		s(5,u_indexode(2));
		s(5,p_indexode(0));

		s(6,x_indexode(4));
		s(6,x_indexode(6));
		s(6,u_indexode(0));
		s(6,u_indexode(1));
		s(6,u_indexode(3));
		s(6,p_indexode(0));

		s(7,x_indexode(7));
		s(7,u_indexode(3));
		s(7,p_indexode(0));

		return true;
	}

	bool ode_diff(DiffStructure& s, double t, const double* x, const double* u, const double* p) override {
	
		return false;
	}

	bool ode_diff_p(DiffStructure& s, double t, const double *x, const double *u, const double *p, int index) override {
	
		return false;

		
		
		return true;
	}


	void u_boundary(double *u_low, double *u_upp) override {

		double umax = 0.1*uscale;

		u_low[0] = -umax;
		u_upp[0] =  umax;
		u_low[1] = -umax;
		u_upp[1] =  umax;
		u_low[2] = -umax;
		u_upp[2] =  umax;
	}

	void x_boundary(double *x_low, double *x_upp) override {

		//x_low[0] = 0;
		//x_upp[0] = 100;

		//x_low[1] = -3;
		//x_upp[1] = 3;

		//x_low[2] = 0; // Additional!!!
		//x_upp[2] = 10;

		//x_low[3] = -10;
		//x_upp[3] = 10;

		//x_low[4] = -4;
		//x_upp[4] = 4;

		//x_low[5] = -1.8325957145940461;
		//x_upp[5] = 1.8325957145940461;

		///x_upp[6] = 102.9;

		//x_low[7] = -10;
		//x_upp[7] = 10;

		//x_low[8] = 0;

	}

	void p_boundary(double *p_low, double *p_upp) override {
		p_low[0] = 1;
		p_upp[0] = 35000;
	}

	void var_boundary(double *x_low, double *x_upp) override {

		double ziel[] =  {6000,200,5,0,0,0,80,0};

		for (int i=0;i<7;i++) {
			x_low[x_index(0,i)] = start[i];
			x_upp[x_index(0,i)] = start[i];
		}
		/*for (int i=0;i<7;i++) {
			x_low[x_index(n_dis-1,i)] = ziel[i];
			x_upp[x_index(n_dis-1,i)] = ziel[i];
		}*/

	}


	void neben_boundary(double *c_low, double *c_upp) override {

		c_low[0] = 6000;
		c_upp[0] = 1e10;

	}

	void neben(double *c, double t, const double *x, const double *u, const double *p) override {

		c[0] = 6000; //sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
	//	double dx[8];
	//	Model3Dof(x, u, dx, c);

	}

	bool neben_structure(DiffStructure &s) override {

		s(0,x_indexode(0));
			s(0,x_indexode(1));
				s(0,x_indexode(2));
	//	s(0,x_index(0,6));
	//	s(0,p_indexode(0));

		return true;
	}
	/*
	bool neben_diff(DiffStructure &s, double t, const double *x, const double *u, const double *p) {

		s(0,x_index(0,4)) =1 ;
		s(0,x_index(0,5)) =-1 ;
		//s.use(0,p_indexode(0));

		return true;
	}

	bool neben_diff_p(DiffStructure &s, double t, const double *x, const double *u, const double *p, int index) {

		s(0,p_indexode(0)) = .01;
		return true;
	}*/

};


/////////////////////////////////////////////////////////////////////////////



int main(int argv, char* argc[]) {

	double mu = 398600.4329;

	double a = 6378+400;
	double epsilon = 0.002;
	double ii = 0.01;
	double Omega = 1;
	double omega = 0.1;
	double phi = 0.;

	double E = atan(sqrt((1-epsilon)/(1+epsilon))*tan(phi/2))*2;
	double M0 = E-epsilon*sin(E);
	double rr = a* (1-epsilon*epsilon)/(1+epsilon*cos(phi));

	Vektor<double> pos(rr*(cos(phi+omega)*cos(Omega)-sin(phi+omega)*cos(ii)*sin(Omega)),
	 rr*(cos(phi+omega)*sin(Omega)+sin(phi+omega)*cos(ii)*cos(Omega)),
	 rr*(sin(phi+omega)*sin(ii)));

	double vv = sqrt(2*mu/rr-mu/a);
	double gamma = atan(epsilon *sin(phi)/(1+epsilon*cos(phi)));
	Vektor<double> vel(vv*(-sin(phi+omega-gamma)*cos(Omega)-cos(phi+omega-gamma)*cos(ii)*sin(Omega)),
		vv*(-sin(phi+omega-gamma)*sin(Omega)+cos(phi+omega-gamma)*cos(ii)*cos(Omega)),
		vv*(cos(phi+omega-gamma)*sin(ii)));


	start[0] = pos.X();
	start[1] = pos.Y();
	start[2] = pos.Z();
	start[3] = vel.X();
	start[4] = vel.Y();
	start[5] = vel.Z();


	TWparameter twparameter("transworhp.xml");
	twparameter.NDIS = 11;
	std::map<std::string,std::string> A = twparameter.Arguments(argv,argc);
	int B =ToInt(A["div"]);
	if (B==0) B=2;

	twparameter.twdiscretization = TWdiscretization(TW_MultipleShooting,0,0);
	twparameter.USERHM=-1;
	twparameter.USERDF=-1;
	twparameter.USERDG=-1;
	
	TWfolder folder(&twparameter,0);

	//XMLNode *xml_laufkatze = TWparameter::ReadParams("flugzeug.xml");
	vector<int> multinode;
	
	for (int i=0;i<B;i++) {
		int a = (double)(i * (twparameter.NDIS-1)) / (B-1);
		multinode.push_back(a);
	}

	OrbitPhase ph( twparameter.NDIS, multinode);

	ph.butcherInit(twparameter.butchertableau,twparameter.stepsize);

	folder.Add(&ph);

	pph=&ph;

	Viewer *viewer = nullptr;
	if (twparameter.PLOT) viewer = new Viewer(&twparameter);

	folder.Init();
	folder.Init(viewer);
	
	//ph.Integrate(6);
	// Ende Optional

	folder.Loop(2);

	delete viewer;

	return 0;
}

