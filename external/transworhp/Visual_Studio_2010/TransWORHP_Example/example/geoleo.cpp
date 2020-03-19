/*----------------------------------------------------------------
 *
 * Geo-Leo-Transfer
 *
 *----------------------------------------------------------------*/

#ifdef WIN32
#include "windows.h"
#endif

#include "TransWORHP.h"

using namespace std;

double geoleoplot (int &len, int &ndgl, int &nsteuer, double *t,double *x, double *u, int &i, int &index);


class GeoLeoPhase : public TransWorhp {
public:

	double rmu ;
	double c;
	double h0;
	double hf;



	GeoLeoPhase(int dis, int START) : TransWorhp("Geo-Leo-Orbit", dis,4,1,1,0,0,0,1) {
		freetime = 1;
		rmu = 62.5;
		c = 1e-2;

		h0 = START;
		hf = 6.6;

		TimeAxis(2.0);



	}
void zen_init(double *zen) {

		// Zen: Referenzwerte der nichtlin. Störungen
		zen[0] = 0;
	}

	void OpenWindows(Viewer *viewer) {

		viewer->PhasePlot("Flight around Earth", geoleoplot,0,1);

	}


	void GetXTitle(int d, char *buf) {
		if (d==0) sprintf(buf,"X_1 Height");
		if (d==1) sprintf(buf,"X_2 Radial Velocity");
		if (d==2) sprintf(buf,"X_3 Tangential Velocity");
		if (d==3) sprintf(buf,"X_4 Polar Angle");

	}

	void GetUTitle(int d, char *buf) {
		if (d==0) sprintf(buf,"U_1 Flight Angle");
	}


	void x_init(double *x, int i, int dis) {



		if (h0==6) {
			x[0] = 6 + .6 * i / (dis-1.);
			x[1] = .03;
			x[2] = .03;
			x[3] = 0 + 9 * i / (dis-1.);
		}
		if (h0<=6) {
			x[0] = 6 + .6 * i / (dis-1.);
			x[1] = .03;
			x[2] = .03;
			x[3] = 0 + 9 * i / (dis-1.);
		}
		if (h0<=2) {
			x[0] = 2 + 5 * i/ (dis-1.);
			x[1] = 0.02;
			x[2] = 5.5 - 2 * i/ (dis-1.);
			x[3] = 0 + 350 * i/ (dis-1.);
		}

		if (h0<=1) {
			x[0] = 1 + 5.5 * i/ (dis-1.);
			x[1] = 0.005;
			x[2] = 8 - 5 * i/ (dis-1.);
			x[3] = 0 + 1400 * i/ (dis-1.);
		}


	}

	void u_init(double *u, int i, int dis) {

		if (h0==6)
			u[0] = -.5 + 1. * i/ (dis-1.);



	}


	void init() {



	}

	void p_init(double *p) {

		if (h0<=6)
			p[0] = 20;

		if (h0<=5)
			p[0] = 40;

		if (h0<=4)
			p[0] = 90;

		if (h0<=3)
			p[0] = 150;

		if (h0<=2)
			p[0] = 250;

		if (h0<=1)
			p[0] = 480;

	}

	double obj() {

		return p(0);
	}


	bool obj_structure(DiffStructure &s) {
		//return false;
		s(0,p_index(0));
		return true;

	}

	bool obj_diff(DiffStructure &s) {
		
		s(0,p_index(0)) = 1; 
	
		return true;
	}

	void ode(double *dx, double t, const double *x, const double *u, const double *p) {

		dx[0] = x[1] * p[0] + ZEN[0];
		dx[1] = (x[2]*x[2]/x[0] - rmu/ (x[0]*x[0]) + c * sin(u[0])) * p[0];
		dx[2] = (-x[1]*x[2]/x[0]                  + c * cos(u[0])) * p[0];
		//!!
		/*double zz = sqrt(1.-u[0]*u[0]);

			dx[1] = (x[2]*x[2]/x[0] - rmu/ (x[0]*x[0]) + c * u[0]) * p[0];
			dx[2] = (-x[1]*x[2]/x[0]                  + c * zz) * p[0];
		*/

		dx[3] = x[2]/x[0]  * p[0];

	}
	bool ode_structure(DiffStructure &s) {


		s(0,x_index(0,1));

		s(1,x_index(0,0));
		s(1,x_index(0,2));
		s(1,u_index(0,0));


		s(2,x_index(0,0));
		s(2,x_index(0,1));
		s(2,x_index(0,2));
		s(2,u_index(0,0));

		s(3,x_index(0,0));
		s(3,x_index(0,2));


		s(0,p_indexode(0));       // p_index-Variante
		s(1,p_indexode(0));       // p_index-Variante
		s(2,p_indexode(0));       // p_index-Variante
		s(3,p_indexode(0));       // p_index-Variante



		return true;
	}


	
bool ode_diff(DiffStructure &s, double t, const double *x, const double *u, const double *p) {
	
		//	dx[0] = x[1] * p[0];
		s(0,x_index(0,1))  = p[0];

		//	dx[1] = (x[2]*x[2]/x[0] - rmu/ (x[0]*x[0]) + c * sin(u[0])) * p[0];
		s(1,x_index(0,0)) = (-x[2]*x[2]/(x[0]*x[0]) + 2.* rmu / (x[0]*x[0]*x[0]) ) * p[0];
		s(1,x_index(0,2)) = 2*x[2]/x[0] * p[0];
		s(1,u_index(0,0)) = c * cos(u[0]) * p[0];

		//	dx[2] = (-x[1]*x[2]/x[0]                  + c * cos(u[0])) * p[0];
		s(2,x_index(0,0)) = x[1]*x[2]/(x[0]*x[0])* p[0] ;
		s(2,x_index(0,1)) = - x[2]/x[0]*p[0]; 
		s(2,x_index(0,2)) = - x[1]/x[0]*p[0];
		s(2,u_index(0,0)) = - c * sin(u[0]) * p[0];

		//	dx[3] = x[2]/x[0]  * p[0];
		s(3,x_index(0,0)) = -x[2]/(x[0]*x[0])*p[0];
		s(3,x_index(0,2)) = 1./x[0]*p[0];

		return true;
	}
	
	bool ode_diff_p(DiffStructure &s, double t, const double *x, const double *u, const double *p, int index) {

		s(0,p_indexode(0)) = x[1] ;
		s(1,p_indexode(0)) = x[2]*x[2]/x[0] - rmu/ (x[0]*x[0]) + c * sin(u[0]);
		s(2,p_indexode(0)) = -x[1]*x[2]/x[0]                  + c * cos(u[0]);
		s(3,p_indexode(0))  = x[2]/x[0];
		
		
		
	
		return true;
	}



	void u_boundary(double *u_low, double *u_upp) {

		u_low[0] = -.2;
		u_upp[0] = .2;

	}

	void x_boundary(double *x_low, double *x_upp) {

		x_low[0] = h0*.9;
		x_upp[0] = hf*1.1;

	}

	void var_boundary(double *x_low, double *x_upp) {

		x_low[x_index(0,0)] = x_upp[x_index(0,0)] = h0;
		x_low[x_index(0,1)] = x_upp[x_index(0,1)] = 0;
		x_low[x_index(0,2)] = x_upp[x_index(0,2)] = sqrt(rmu/h0);
		x_low[x_index(0,3)] = x_upp[x_index(0,3)] = 0;

		x_low[x_index(n_dis-1,0)] = x_upp[x_index(n_dis-1,0)] = hf;         // GEO
		x_low[x_index(n_dis-1,1)] = x_upp[x_index(n_dis-1,1)] = 0;
		x_low[x_index(n_dis-1,2)] = x_upp[x_index(n_dis-1,2)] = sqrt(rmu/hf);

	}

	void p_boundary(double *p_low, double *p_upp) {

		p_low[0] = 1;

		if (h0<=2) p_low[0] = 150;
		if (h0<=1) p_low[0] = 400;

		p_upp[0] = 2000;
	}


};




/*	void terminate() {

		cout << "t_f      = " << p(0) << endl;
		cout << "x_4(t_f) = " << x(n_dis-1,3) << " = " << x(n_dis-1,3) /M_PI << "pi" << endl;
		
		cout << worhp_c->status << endl;
	
	}*/


double geoleoplot (int &len, int &ndgl, int &nsteuer, double *t,double *x, double *u, int &i, int &index) {

	//cout << x[p0->x_index(i,0)] << " " << x[p0->x_index(i,1)] << endl;
	if (index==0) return x[i*(ndgl+nsteuer)] * sin(x[i*(ndgl+nsteuer)+3]);
	if (index==1) return x[i*(ndgl+nsteuer)] * cos(x[i*(ndgl+nsteuer)+3]);

}


/////////////////////////////////////////////////////////////////////////////



int main(int argv, char* argc[]) {
/*
	TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv,argc);
	TWfolder folder(&twparameter,0);

	double START = 4;

	twparameter.NDIS=51;

	GeoLeoPhase ph(twparameter.NDIS,START);
	folder.Add(&ph);

	Viewer *viewer = 0;
	if (twparameter.PLOT) viewer = new Viewer(&twparameter);

	folder.Init();
	folder.Init(viewer);

	ph.Integrate();

	folder.Loop();
	
	
	cout << "t_f      = " << ph.p(0) << endl;
		cout << "x_4(t_f) = " << ph.x(ph.n_dis-1,3) << " = " << ph.x(ph.n_dis-1,3) /M_PI << "pi" << endl;
		
		cout << folder.worhp_c.status << endl;

	delete viewer;
	return 0;
	*/
	
	
	
	
	TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv,argc);

	double START = 3;

	Viewer *viewer = 0;
	if (twparameter.PLOT) viewer = new Viewer(&twparameter);

	// Startschaetzung in Datei erstellen
	{
		twparameter.NDIS=51;
		GeoLeoPhase ph(twparameter.NDIS,START);
		TWfolder folder(&twparameter,0);
		folder.Add(&ph);
		folder.Init();
		folder.Init(viewer);

		folder.Loop(1);
		ph.ToMATLAB("geoleo51.m");
	}

	if (viewer) viewer->CloseAll();

	// Startschaetzung aus Datei lesen (und fuer mehr Diskretisierungspunkte interpolieren)
	{
		twparameter.NDIS=20001;
		GeoLeoPhase ph(twparameter.NDIS,START);
		TWfolder folder(&twparameter,0);
		folder.Add(&ph);
		folder.Init();
		folder.Init(viewer);

		ph.FromMATLAB("geoleo51.m");
		folder.Loop(1);
		
		
		
		
		
		
		
	/*		cout << "Writing matrix files..." << endl;
	TransWorhp::DoubleToMATLAB(folder.worhp_o.X, folder.worhp_o.n, "X_GEOLEO.m");
	TransWorhp::DoubleToMATLAB(ph.T, ph.n_dis, "T_GEOLEO.m");
	
	char var[] = {'X','G','F','L','M',0};
	char pert[] = {'P','R','Q','B',0};
	int i=0;
	while (var[i]) {
		//os << "DX/DP" << endl;
		ofstream os( (string("ZenD") + var[i] + ".dat").c_str());
		
		int j=0;
		while (pert[j]) {
			//os << "DX/DP" << endl;
			folder.WriteSensitivity(os, var[i], pert[j]);
			j++;
		}
		
		i++;
	}
*/
	}







	delete viewer;
	
	return 0;
	
}

