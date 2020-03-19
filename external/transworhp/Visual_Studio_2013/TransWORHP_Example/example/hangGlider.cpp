/*----------------------------------------------------------------
 *
 * Hang Glider aus Betts S.283ff
 * maximiere Reichweite
 *
 *----------------------------------------------------------------*/

#ifdef WIN32
#include "windows.h"
#endif

#include "TransWORHP.h"
#include "../src/core/TWparameter.h"
#include "../src/core/Viewer.h"

using namespace std;


class HangGliderPhase : public TransWorhp {
public:

	double uM, R, C0, k, m, S, rho, g;
	HangGliderPhase(int dis) : TransWorhp("Hang Glider",dis,4,1,1,0,0) {
	
	freetime = 1;
	  
	uM = 2.5;
	R = 100.0;
	C0 = 0.034;
	k = 0.069662;
	m = 100.0;
	S = 14.0;
	rho = 1.13;
	g = 9.80665;
	
	}
	

	void x_init(double *x, int i, int dis) {
		/*
		x[0] = 0;
		x[1] = 1000;
		x[2] = 13.23;
		x[3] = -1.29;
		*/
		
		x[0] = 1000-1000*(dis-i)/dis;
		x[1] = 100*(dis-i)/dis+900;
		x[2] = 13.23;
		x[3] = -1.29;
	}

	void u_init(double *u, int i, int dis) {

		u[0] = 1.0;
	}
	
	void p_init(double *p) {
		p[0] = 100.0;
	}

	double obj() {
		return -x(n_dis-1, 0);
	}


	bool obj_structure(DiffStructure &s) {

		s(0,x_index(n_dis-1,0));
		return true;
	}

	bool obj_diff(DiffStructure &s) {

		s(0,x_index(n_dis-1,0)) = -1.0;
		return false;
	}
	
	void ode(double *dx, double t, const double *x, const double *u, const double *p) {

		double CD;
		double D, L, X, ua, Vy, vr, sinn, cosn;

		CD = C0+k*u[0]*u[0];
		
		
		X = (x[0]/R-2.5)*(x[0]/R-2.5);
		
		ua = uM*(1-X)*exp(-X);
		
		Vy = x[3] - ua;
		
		vr = sqrt(x[2]*x[2] + Vy*Vy);
		
		sinn = Vy/vr;
		cosn = x[2]/vr;
		
		D = 0.5*CD*rho*S*vr*vr;
		
		L = 0.5*u[0]*rho*S*vr*vr;
		
	  
		dx[0] = x[2];
		dx[1] = x[3];
		dx[2] = 1/m*( -L * sinn -D * cosn );
		dx[3] = 1/m*( L * cosn - D * sinn -m*g );
		
		dx[0] *= p[0];
		dx[1] *= p[0];
		dx[2] *= p[0];
		dx[3] *= p[0];
		
		//cout << dx[0] << " " << dx[1] << " " << dx[2] << " " << dx[3] << endl;
	}

	bool ode_structure(DiffStructure &s) {
		
		
		return false;
	}

	bool ode_diff(DiffStructure &s, double t, const double *x, const double *u, const double *p) {
		
		
		return false;
	}
	
	bool ode_diff_p(DiffStructure& s, double t, const double *x, const double *u, const double *p, int index) {
	
		return false;
	}

	void u_boundary(double *u_low, double *u_upp) {
		u_low[0] = 0.0;
		u_upp[0] = 1.4;
	}

	void x_boundary(double *x_low, double *x_upp) {

		x_low[0] = 0.0;

	}
	
	void p_boundary(double *p_low, double *p_upp) {
		p_low[0] = 80.0;
		p_upp[0] = 120.0;
	}

	void var_boundary(double *x_low, double *x_upp) {

		x_low[x_index(0,0)] = x_upp[x_index(0,0)] = 0.0;
		x_low[x_index(0,1)] = x_upp[x_index(0,1)] = 1000.0;
		x_low[x_index(0,2)] = x_upp[x_index(0,2)]= 13.227567500;
		x_low[x_index(0,3)] = x_upp[x_index(0,3)] = -1.2875005200;
		
		//x_low[x_index(n_dis-1,0)] = x_upp[x_index(n_dis-1,0)] = ; //frei
		x_low[x_index(n_dis-1,1)] = x_upp[x_index(n_dis-1,1)] = 900.0;
		x_low[x_index(n_dis-1,2)] = x_upp[x_index(n_dis-1,2)] = 13.227567500;
		x_low[x_index(n_dis-1,3)] = x_upp[x_index(n_dis-1,3)] = -1.2875005200;
	}
	
	void terminate() {

		vector<double> schrittweite;
		for (int i = 0; i < n_dis; i++) {
			schrittweite.push_back(T[i]);
		}
		SCHRITTWEITE.push_back(schrittweite);
		schrittweite.clear();

		diskretisierungsfehler(); // <<---- Aufruf fuer Anpassung
	}
	
	void OpenWindows(Viewer *viewer) {
		viewer->disFehler(SCHRITTWEITE, FEHLER,twdiscretization);
		//viewer->disFehler(SCHRITTWEITE, FEHLER2);
		viewer->restrictionPlot(T,Lambda,n_dis,n_ctrl,n_ode,twdiscretization);
		viewer->adjungiertenPlot(T,Mu,n_dis,n_ctrl,n_ode,n_con,twdiscretization);
	}
	
};

/////////////////////////////////////////////////////////////////////////////

int main(int argv, char* argc[]) {

	TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv, argc);

	vector<vector<double> > schrittweite;
	vector<vector<double> > fehler;
	vector<double> T_alt;
	vector<double> T_neu;

	int alt;
	int neu;


	const double tol = twparameter.meshref_tol;
	double max = 1.0;

	Viewer *viewer = 0;
	

	if (twparameter.PLOT) viewer = new Viewer(&twparameter);


	// Startschaetzung in Datei erstellen
	{
		HangGliderPhase ph(twparameter.NDIS);

		TWfolder folder(&twparameter, 0);
		folder.Add(&ph);
		folder.Init();
		folder.Init(viewer); //<------------------------- grafik ausgabe
		
		
		folder.Loop(1);
		
		ph.ToMATLAB("hangGliderTEST.m");
		
		T_alt = ph.getGrid();
		T_neu = ph.refineAlg(max);
		
		neu = T_neu.size();
		
		schrittweite = ph.SCHRITTWEITE;
		fehler = ph.FEHLER;
		
	}
	

	int k = 0;

	while (max > tol) {
	  
		/*if (k%2 == 0) {
			twparameter.twdiscretization = TWdiscretization(TW_HermiteSimpson,2,1);
		} else {
			twparameter.twdiscretization = TWdiscretization(TW_Trapez,1,0);
		}*/    

		if (viewer) viewer->CloseAll();

		HangGliderPhase ph(neu);
		
		ph.newGrid(T_neu);

		ph.SCHRITTWEITE = schrittweite;
		ph.FEHLER = fehler;

		TWfolder folder(&twparameter, 0);
		folder.Add(&ph);
		folder.Init();

		
		folder.Init(viewer);
		
		ph.FromMATLAB("hangGliderTEST.m");

		folder.Loop(1);

		ph.ToMATLAB("hangGliderTEST.m");

		schrittweite = ph.SCHRITTWEITE;
		fehler = ph.FEHLER;

		alt = neu;

		T_alt = T_neu;
		T_neu = ph.refineAlg(max);
		neu = T_neu.size();


		k++;
	}

	if (viewer) viewer->CloseAll();

	cout << "........ fertig nach " << k << " Schritten ........" << endl;

	// Endwert anzeigen
	{
		HangGliderPhase ph(neu);

		ph.newGrid(T_neu);

		ph.SCHRITTWEITE = schrittweite;
		ph.FEHLER = fehler;

		TWfolder folder(&twparameter, 0);
		folder.Add(&ph);
		folder.Init();

		
		folder.Init(viewer);
		
		ph.FromMATLAB("hangGliderTEST.m");

		folder.Loop(0,1);

		if (viewer) viewer->CloseAll();

	}
	
	delete viewer;

	return 0;
}

