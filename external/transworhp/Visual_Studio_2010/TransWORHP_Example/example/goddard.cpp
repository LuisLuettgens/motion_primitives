/*----------------------------------------------------------------
 *
 *  Example: Goddard Problem / Hoehenrakete
 *
 *----------------------------------------------------------------*/

#ifdef WIN32
#include "windows.h"
#endif

#include "../src/core/TransWORHP.h"

using namespace std;


class Goddard : public TransWorhp {
public:

	Goddard(int dis) : TransWorhp("Hoehenrakete", dis, 3, 1, 1, 0, 0) {
		freetime = 1;
	}

	string GetXTitle(int d) {
		if (d==0) return "h - Hoehe";
		if (d==1) return "v - Geschwindigkeit";
		if (d==2) return "m - Masse";
	}

	string GetUTitle(int d) {
		if (d==0) return "u - Steuerung";
	}

        void x_init(double *x, int i, int dis) {
		x[0] = 1.0;
		x[1] = 1.0;
		x[2] = 1.0;
	}

	void p_init(double *p) {
		p[0] = 211.6;
	}

	void u_init(double *u, int i, int dis) {
		u[0] = 1.0;
	}

	double obj() {
		return -x(n_dis-1, 0);
	}

	bool obj_structure(DiffStructure &s) {
		s(0, x_index(n_dis-1, 0));
		return true;
	}

	bool obj_diff(DiffStructure &s) {
		s(0, x_index(n_dis-1, 0) ) = -1.0;
		return true;
	}

	void ode(double *dx, double t, const double *x, const double *u, const double *p) {
		
		const double c = 2060.103;
		const double alpha = 1.22710e-2;
		const double beta = 1.45010e-4;
		const double g0 = 9.810;
		const double r0 = 6.37110e6;
		
		double g = g0 * (r0/(r0+x[0]))*(r0/(r0+x[0]));
		double d = alpha*x[1]*x[1] * exp(-beta*x[0]);
	  
		dx[0] = x[1];
		dx[1] = 1.0/x[2]*(u[0]*c - d) - g;
		dx[2] = -u[0];

		// free final time
		dx[0] *= p[0];
		dx[1] *= p[0];
		dx[2] *= p[0];
	}

	bool ode_structure(DiffStructure &s) {
		/*
		s(0, x_indexode(1));

		s(1, x_indexode(0));
		s(1, x_indexode(1));
		s(1, x_indexode(2));
		s(1, u_indexode(0));

		s(2, u_indexode(0));
		
		s(0, p_indexode(0));
		s(1, p_indexode(0));
		s(2, p_indexode(0));
		*/
		return false;
	}

	bool ode_diff(DiffStructure& s, double t, const double* x, const double* u, const double* p) {
		
		/*
		s(0, x_indexode(1)) = 1.0;

		s(1, x_indexode(0));
		s(1, x_indexode(1));
		s(1, x_indexode(2));
		s(1, u_indexode(0));

		s(2, u_indexode(0)) = -1.0;
		
		s(0, p_indexode(0)) = x[1];
		s(1, p_indexode(0)) = 1.0/x[2]*(u[0]*c - d) - g;
		s(2, p_indexode(0)) = -u[0];
		*/
		
		return false;
	}

	bool ode_diff_p(DiffStructure& s, double t, const double *x, const double *u, const double *p, int index) {
		
		return false;
	}

	void x_boundary(double *x_low, double *x_upp) {
	  
	}

	void u_boundary(double *u_low, double *u_upp) {
		u_low[0] = 0.0;
		u_upp[0] = 10.0;
	}

	void p_boundary(double *p_low, double *p_upp) {
		p_low[0] = 200.0;
		p_upp[0] = 300.0;
	}

	void var_boundary(double *x_low, double *x_upp) {

		x_low[x_index(0,0)] = x_upp[x_index(0,0)] = 0.0;
		x_low[x_index(0,1)] = x_upp[x_index(0,1)] = 0.0;
		x_low[x_index(0,2)] = x_upp[x_index(0,2)] = 214.839;

		x_low[x_index(n_dis-1,2)] = x_upp[x_index(n_dis-1,2)] = 67.9833;
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
		Goddard ph(twparameter.NDIS);

		TWfolder folder(&twparameter, 0);
		folder.Add(&ph);
		folder.Init();
		folder.Init(viewer); //<------------------------- grafik ausgabe
		
		
		folder.Loop(1);
		
		ph.ToMATLAB("goddardTEST.m");
		
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

		Goddard ph(neu);
		
		ph.newGrid(T_neu);

		ph.SCHRITTWEITE = schrittweite;
		ph.FEHLER = fehler;

		TWfolder folder(&twparameter, 0);
		folder.Add(&ph);
		folder.Init();

		
		folder.Init(viewer);
		

		ph.FromMATLAB("goddard.m");

		folder.Loop(1);

		ph.ToMATLAB("goddardTEST.m");

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
		Goddard ph(neu);

		ph.newGrid(T_neu);

		ph.SCHRITTWEITE = schrittweite;
		ph.FEHLER = fehler;

		TWfolder folder(&twparameter, 0);
		folder.Add(&ph);
		folder.Init();

		
		folder.Init(viewer);
		
		ph.FromMATLAB("goddardTEST.m");

		folder.Loop(0,1);

		if (viewer) viewer->CloseAll();

	}
	
	delete viewer;

	return 0;
}
