/*----------------------------------------------------------------
 *
 *  Example: Splineproblem
 *
 *----------------------------------------------------------------*/

#ifdef WIN32
#include "windows.h"
#endif

#include "TransWORHP.h"

#include "../src/core/TWparameter.h"
#include "../src/core/Viewer.h"

using namespace std;

class SplinePhase : public TransWorhp {
public:

	double alpha;
	
	SplinePhase(int dis, double alpha) : TransWorhp("Spline",dis,3,1,0,0,0), alpha(alpha) {}

	string GetXTitle(int d) {
		if (d==0) return "x_1";
		if (d==1) return "x_2";
		if (d==2) return "x_3";
	}

	string GetUTitle(int d) {
		if (d==0) return "u";
	}
	

	double obj() {
		return .5 * x(n_dis-1,2);
	}

	bool obj_structure(DiffStructure &s) {
		s(0, x_index(n_dis-1,2));
		return true;
	}

	bool obj_diff(DiffStructure &s) {
		s(0,x_index(n_dis-1,2)) = .5;
		return true;
	}


	void ode(double *dx, double t, const double *x, const double *u, const double *p) {
		dx[0] = x[1];
		dx[1] = u[0];
		dx[2] = u[0]*u[0];
	}

	bool ode_structure(DiffStructure &s) {
		s(0, x_indexode(1)); // dx[0] / dx[1]
		s(1, u_indexode(0)); // dx[1] / du[0]
		s(2, u_indexode(0)); // dx[2] / du[0]
		return true;
	}

	bool ode_diff(DiffStructure &s, double t, const double *x, const double *u, const double *p) {
		s(0, x_indexode(1))= 1;
		s(1, u_indexode(0))= 1;
		s(2, u_indexode(0))= 2*u[0];
		return true;
	}


	void x_boundary(double *x_low, double *x_upp) {
		x_upp[0] = alpha;
	}
	
	void u_boundary(double *u_low, double *u_upp) {
		//u_upp[0] = 0;
	}

	void var_boundary(double *x_low, double *x_upp) {
		x_low[x_index ( 0,0 ) ] = x_upp[x_index ( 0,0 ) ] = 0;
		x_low[x_index ( 0,1 ) ] = x_upp[x_index ( 0,1 ) ] = 1;
		x_low[x_index ( 0,2 ) ] = x_upp[x_index ( 0,2 ) ] = 0;

		x_low[x_index ( n_dis-1,0 ) ] = x_upp[x_index ( n_dis-1,0 ) ] = 0;
		x_low[x_index ( n_dis-1,1 ) ] = x_upp[x_index ( n_dis-1,1 ) ] = -1;
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
		viewer->disFehler(SCHRITTWEITE, FEHLER, twdiscretization);
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

	
	// Parameter alpha aus Kommandozeile lesen
	double alpha = 1/9.;
	
	
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
		SplinePhase ph(alt, alpha);

		TWfolder folder(&twparameter, 0);
		folder.Add(&ph);
		folder.Init();
		folder.Init(viewer); //<------------------------- grafik ausgabe

		folder.Loop(1);
		

		ph.ToMATLAB("splineTEST.m");
		ph.ToMATLAB_LambdaMu("splineLM.m");
		
		T_alt = ph.getGrid();
		T_neu = ph.refineAlg(max);

		neu = T_neu.size();

		schrittweite = ph.SCHRITTWEITE;
		fehler = ph.FEHLER;
	}
	


	int k = 0;

	while (max > tol) {

		if (viewer) viewer->CloseAll();


		SplinePhase ph(neu, alpha);

		ph.newGrid(T_neu);

		ph.SCHRITTWEITE = schrittweite;
		ph.FEHLER = fehler;

		TWfolder folder(&twparameter, 0);
		folder.Add(&ph);
		folder.Init();

		
		//folder.Init(viewer);
		
		
		ph.FromMATLAB("splineTEST.m");
		ph.FromMATLAB_LambdaMu("splineLM.m");

		folder.Loop(1,0);

		ph.ToMATLAB("splineTEST.m");
		ph.ToMATLAB_LambdaMu("splineLM.m");

		schrittweite = ph.SCHRITTWEITE;
		fehler = ph.FEHLER;

		alt = neu;

		T_alt = T_neu;

		
		T_neu = ph.refineAlg(max);


		neu = T_neu.size();
		

		k++;

	}
	

	cout << "........ fertig nach " << k << " Schritten ........" << endl;

	if (viewer) viewer->CloseAll();

	// Endwert anzeigen
	{
		SplinePhase ph(neu, alpha);

		ph.newGrid(T_neu);

		ph.SCHRITTWEITE = schrittweite;
		ph.FEHLER = fehler;

		TWfolder folder(&twparameter, 0);
		folder.Add(&ph);
		folder.Init();

		
		folder.Init(viewer);
		

		ph.FromMATLAB("splineTEST.m");
		ph.FromMATLAB_LambdaMu("splineLM.m");
		
		folder.Loop(1,0);

	}

	end = clock();

	cout << end - start << "  -  Schritte:" << k + 1 << endl;


	delete viewer;

	return 0;
}

