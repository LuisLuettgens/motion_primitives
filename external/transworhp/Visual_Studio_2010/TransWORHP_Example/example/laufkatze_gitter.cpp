/*----------------------------------------------------------------
 *
 * Laufkatze: Einfache Trajektorie
 *
 *----------------------------------------------------------------*/

#ifdef WIN32
#include "windows.h"
#endif

#include "TransWORHP.h"

#include "TWspline.h"

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


class LaufkatzePhase : public TransWorhp {
public:
	double grav;
	XMLNode* scenenode;

	LaufkatzePhase(XMLNode *xmlmain, int dis) : TransWorhp("Laufkatze", dis,9,2,1,0,0) {
		freetime = 1;
		grav = 9.81;

		if (xmlmain)
			scenenode = xmlmain->GetFirstChild("SCENE");

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
		p[0] = 12;
	}

	void x_init(double *x, int i, int dis) {
		x[5] = 5;
	}
	void u_init(double *u, int i, int dis) {
		u[0] = 0;
		u[1] = 0;
	}


	double obj() {
		return x(n_dis-1,8) + .1 * p(0);
	}


	bool obj_structure(DiffStructure &s) {
		s(0,x_index(n_dis-1,8) );
		s(0,p_index(0) );
		return true;

	}

	bool obj_diff(DiffStructure &s) {
		s(0,x_index(n_dis-1,8)) = 1 ;
		s(0,p_index(0)) = .1 ;
		return true;
	}


	void ode(double *dx, double t, const double *x, const double *u, const double *p) {

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
	bool ode_structure(DiffStructure &s) {

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

		x_low[5] = 0.5;
		x_upp[5] = 15;

		x_low[6] = -3;
		x_upp[6] = 3;

		x_low[7] = -10;
		x_upp[7] = 10;

		x_low[8] = 0;

	}

	void p_boundary(double *p_low, double *p_upp) {
		p_low[0] = 1;
		p_upp[0] = 15;
	}

	void var_boundary(double *x_low, double *x_upp) {

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


	void terminate() {

		vector<double> schrittweite;
		for (int i = 0; i < n_dis; i++) {
			schrittweite.push_back(T[i]);
		}
		SCHRITTWEITE.push_back(schrittweite);
		schrittweite.clear();

		diskretisierungsfehler(); // muss aufgerufen werden, damit Fehler geplottet wird
		
		//diskretisierungsfehler2();
	}

	void OpenWindows(Viewer *viewer) {

		//viewer->ThreeD("Hallenansicht", scenenode, laufkatze3d);

		twparameter;
		viewer->disFehler(SCHRITTWEITE, FEHLER,twdiscretization);
		//viewer->disFehler(SCHRITTWEITE, FEHLER2);
		viewer->restrictionPlot(T,Lambda,n_dis,n_ctrl,n_ode,twdiscretization);
		viewer->adjungiertenPlot(T,Mu,n_dis,n_ctrl,n_ode,n_con,twdiscretization);

	}

};


/////////////////////////////////////////////////////////////////////////////



int main(int argv, char* argc[]) {

	clock_t start, end;

	start = clock();

	TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv, argc);

	XMLNode *xml_laufkatze = TWparameter::ReadParams("laufkatze.xml");

	vector<vector<double> > schrittweite;
	vector<vector<double> > fehler;
	vector<double> T_alt;
	vector<double> T_neu;

	int alt;
	int neu;


	const double tol = twparameter.meshref_tol;
	double max = 1;

	Viewer *viewer = 0;

	if (twparameter.PLOT) viewer = new Viewer(&twparameter);


	// Startschaetzung in Datei erstellen
	{
		LaufkatzePhase ph(xml_laufkatze, twparameter.NDIS);


		TWfolder folder(&twparameter, 0);
		folder.Add(&ph);
		folder.Init();
		folder.Init(viewer); //<------------------------- grafik ausgabe

		folder.meshRef();
		
		/*
		folder.Loop(0);

		ph.ToMATLAB("katzeTEST.m");

		T_alt = ph.getGrid();
		T_neu = ph.refineAlg(max, twparameter);

		neu = T_neu.size();

		schrittweite = ph.SCHRITTWEITE;
		fehler = ph.FEHLER;
		*/
	}

	/*

	int k = 0;

	while (max > tol) {
	  
		if (k > 2) {
			twparameter.twdiscretization = TWdiscretization(TW_HermiteSimpson,2,1);
		} else {
			twparameter.twdiscretization = TWdiscretization(TW_Trapez,1,0);
		}

		if (viewer) viewer->CloseAll();


		LaufkatzePhase ph(xml_laufkatze, neu);

		ph.newGrid(T_neu);

		ph.SCHRITTWEITE = schrittweite;
		ph.FEHLER = fehler;

		TWfolder folder(&twparameter, 0);
		folder.Add(&ph);
		folder.Init();

		
		//folder.Init(viewer);
		

		ph.FromMATLAB("katzeTEST.m");

		folder.Loop(0);

		ph.ToMATLAB("katzeTEST.m");

		schrittweite = ph.SCHRITTWEITE;
		fehler = ph.FEHLER;


		alt = neu;

		T_alt = T_neu;
		T_neu = ph.refineAlg(max, twparameter);
		neu = T_neu.size();
  
		
		cerr << "max: " << max << endl;

		k++;
	}



	cout << "........ fertig nach " << k << " Schritten ........" << endl;

	// Endwert anzeigen
	{
		if (viewer) viewer->CloseAll();
	  
		LaufkatzePhase ph(xml_laufkatze, neu);

		ph.newGrid(T_neu);

		ph.SCHRITTWEITE = schrittweite;
		ph.FEHLER = fehler;

		TWfolder folder(&twparameter, 0);
		folder.Add(&ph);
		folder.Init();

		
		folder.Init(viewer);
		

		ph.FromMATLAB("katzeTEST.m");

		folder.Loop(0, 0);

		if (viewer) viewer->CloseAll();


		end = clock();

		cout << end - start << "  -  Schritte:" << k + 1 << endl;
	}

	int a;
	std::cin >> a;

	*/

	delete viewer;


	/*
	TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv,argc);
	TWfolder folder(&twparameter,0);

	XMLNode *xml_laufkatze = TWparameter::ReadParams("laufkatze.xml");



	LaufkatzePhase ph(xml_laufkatze, twparameter.NDIS);
	folder.Add(&ph);

	Viewer *viewer = 0;
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
	
	
	ph.ToMATLAB("referenz.dat");
	delete viewer;
	delete xml_laufkatze;
	*/

	return 0;
}

