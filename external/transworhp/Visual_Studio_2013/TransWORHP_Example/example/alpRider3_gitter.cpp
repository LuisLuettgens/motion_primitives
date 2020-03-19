/*----------------------------------------------------------------
*
* Alp Rider aus Betts
*
*----------------------------------------------------------------*/

#ifdef WIN32
#include "windows.h"
#endif

// Zeitmessung
#include <time.h>

#include "TransWORHP.h"

using namespace std;


class AlpRiderPhase : public TransWorhp {
public:

	AlpRiderPhase(int dis) : TransWorhp("Alp Rider", dis, 4, 2, 0, 0, 1, 1) {
		//	freetime = 1;
		lagrange_weight[0] = 1;
	}


	double peak(double t, double a, double b) {
		return exp(-b*(t - a)*(t - a));
	}

	/*
	void GetXTitle(int d, char *s) {

	if (d==0) strcpy(s, "Zustand 1");
	if (d==1) strcpy(s, "Zustand 2");
	if (d==2) strcpy(s, "Zustand 3");
	if (d==3) strcpy(s, "Zustand 4");
	if (d==4) strcpy(s, "Integral");

	}

	void GetUTitle(int d, char *s) {

	if (d==0) strcpy(s, "Steuerung 1");
	if (d==1) strcpy(s, "Steuerung 2");

	}
	*/


	void x_init(double *x, int i, int dis) {
		x[0] = 2;
		x[1] = 1;
		x[2] = 1;
		x[3] = 1;
	}
	void u_init(double *u, int i, int dis) {
		u[0] = 1;
		u[1] = 1;
	}


	double obj() {
		return 0;//x(n_dis-1,4);
	}


	bool obj_structure(DiffStructure &s) {
		//	s(0,x_index(n_dis-1,4) );
		return true;

	}
	/*
	bool obj_diff(DiffStructure &s) {
	s(0,x_index(n_dis-1,8)) = 1 ;
	s(0,p_index(0)) = .1 ;
	return true;
	}
	*/


	void integral(double *f, double t, const double *x, const double *u,
		const double *p) {

		double c = 20;

		f[0] = ((100.0 * (x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3]) +
			0.01 * (u[0] * u[0] + u[1] * u[1]))) * c;

	}

	bool integral_structure(DiffStructure &s) {

		s(0, x_indexode(0));
		s(0, x_indexode(1));
		s(0, x_indexode(2));
		s(0, x_indexode(3));
		s(0, u_indexode(0));
		s(0, u_indexode(1));
		return true;
	}

	bool integral_diff(DiffStructure &s, double t, const double *x, const double *u,
		const double *p) {

		double c = 20;


		s(0, x_indexode(0)) = c*200.0 * x[0];
		s(0, x_indexode(1)) = c*200.0 * x[1];
		s(0, x_indexode(2)) = c*200.0 * x[2];
		s(0, x_indexode(3)) = c*200.0 * x[3];
		s(0, u_indexode(0)) = c*0.02 * u[0];
		s(0, u_indexode(1)) = c*0.02 * u[1];


		return true;
	}


	void ode(double *dx, double t, const double *x, const double *u, const double *p) {

		double c = 20;
		dx[0] = (-10.0 * x[0] + u[0] + u[1]) * c;
		dx[1] = (-2.0 * x[1] + u[0] + 2.0 * u[1]) * c;
		dx[2] = (-3.0 * x[2] + 5.0 * x[3] + u[0] - u[1]) * c;
		dx[3] = (5.0 * x[2] - 3.0 * x[3] + u[0] + 3.0 * u[1]) * c;
	}


	// Optional

	bool ode_structure(DiffStructure &s) {

		s(0, x_indexode(0));
		s(0, u_indexode(0));
		s(0, u_indexode(1));

		s(1, x_indexode(1));
		s(1, u_indexode(0));
		s(1, u_indexode(1));

		s(2, x_indexode(2));
		s(2, x_indexode(3));
		s(2, u_indexode(0));
		s(2, u_indexode(1));

		s(3, x_indexode(2));
		s(3, x_indexode(3));
		s(3, u_indexode(0));
		s(3, u_indexode(1));

		return true;
	}


	bool ode_diff(DiffStructure& s, double t, const double* x, const double* u, const double* p) {

		double c = 20;

		//dx[0] = (-10.0 * x[0] + u[0] + u[1]) * p[0];
		s(0, x_indexode(0)) = -10.0 * c;
		s(0, u_indexode(0)) = c;
		s(0, u_indexode(1)) = c;

		//dx[1] = (-2.0 * x[1] + u[0] + 2.0 * u[1]) * p[0];
		s(1, x_indexode(1)) = -2.0 * c;
		s(1, u_indexode(0)) = c;
		s(1, u_indexode(1)) = 2.0 * c;

		//dx[2] = (-3.0 * x[2] + 5.0 * x[3] + u[0] - u[1]) * p[0];
		s(2, x_indexode(2)) = -3.0 * c;
		s(2, x_indexode(3)) = 5.0 * c;
		s(2, u_indexode(0)) = c;
		s(2, u_indexode(1)) = -c;

		//dx[3] = (5.0 * x[2] - 3.0 * x[3] + u[0] + 3.0 * u[1]) * p[0];
		s(3, x_indexode(2)) = 5.0 * c;
		s(3, x_indexode(3)) = -3.0 * c;
		s(3, u_indexode(0)) = c;
		s(3, u_indexode(1)) = 3.0 * c;

		return true;
	}





	void u_boundary(double *u_low, double *u_upp) {

	}

	void x_boundary(double *x_low, double *x_upp) {

	}


	void var_boundary(double *x_low, double *x_upp) {

		x_low[x_index(0, 0)] = 2;
		x_upp[x_index(0, 0)] = 2;
		x_low[x_index(n_dis - 1, 0)] = 2;
		x_upp[x_index(n_dis - 1, 0)] = 2;

		x_low[x_index(0, 1)] = 1;
		x_upp[x_index(0, 1)] = 1;
		x_low[x_index(n_dis - 1, 1)] = 3;
		x_upp[x_index(n_dis - 1, 1)] = 3;

		x_low[x_index(0, 2)] = 2;
		x_upp[x_index(0, 2)] = 2;
		x_low[x_index(n_dis - 1, 2)] = 1;
		x_upp[x_index(n_dis - 1, 2)] = 1;

		x_low[x_index(0, 3)] = 1;
		x_upp[x_index(0, 3)] = 1;
		x_low[x_index(n_dis - 1, 3)] = -2;
		x_upp[x_index(n_dis - 1, 3)] = -2;

	}

	void neben(double *c, double t, const double *x, const double *u, const double *p) {
		double cc = 20;

		c[0] = x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3] - 3.0*peak(t*cc,
			3, 12) - 3.0*peak(t*cc, 6, 10) - 3.0*peak(t*cc, 10, 6) -
			8.0*peak(t*cc, 15, 4);
		//c[0] = x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3] - 3.0*peak(t*p[0], 3, 12);
	}

	void neben_boundary(double *c_low, double *c_upp) {
		c_low[0] = 0.01;
	}

	bool neben_structure(DiffStructure &s) {

		//s(0, x_index(0, 0));
		//s(0, x_index(0, 1));
		//s(0, x_index(0, 2));
		//s(0, x_index(0, 3));
		return false;
	}

	bool neben_diff(DiffStructure &s, double t, const double *x, const double *u, const double *p) {

		// c[0] = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3] - 3.0*peak(t*p[0], 3, 12) - 3.0*peak(t*p[0], 6, 10) - 3.0*peak(t*p[0], 10, 6) - 8.0*peak(t*p[0], 15, 4);
		//s(0, x_index(0, 0)) = 2.0 * x[0];
		//s(0, x_index(0, 1)) = 2.0 * x[1];;
		//s(0, x_index(0, 2)) = 2.0 * x[2];;
		//s(0, x_index(0, 3)) = 2.0 * x[3];;

		return false;
	}


	void terminate() {

		vector<double> schrittweite;
		for (int i = 0; i < n_dis; i++) {
			schrittweite.push_back(T[i]);
		}
		SCHRITTWEITE.push_back(schrittweite);
		schrittweite.clear();

		diskretisierungsfehler(); // muss aufgerufen werden, damit Fehler geplottet wird
	}

	void OpenWindows(Viewer *viewer) {
		viewer->disFehler(SCHRITTWEITE, FEHLER,twdiscretization);
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

	vector<vector<double> > schrittweite;
	vector<vector<double> > fehler;
	vector<double> T_alt;
	vector<double> T_neu;

	int alt = 21; //Anz der Startpunkte
	int neu;


	const double tol = twparameter.meshref_tol;
	double max = 1;

	Viewer *viewer = 0;

	if (twparameter.PLOT) viewer = new Viewer(&twparameter);


	// Startschaetzung in Datei erstellen
	{
		AlpRiderPhase ph(201);


		TWfolder folder(&twparameter, 0);
		folder.Add(&ph);
		folder.Init();
		folder.Init(viewer); //<------------------------- grafik ausgabe


		folder.Loop(0);
/*

		ph.ToMATLAB("alpTEST.m");

		T_alt = ph.getGrid();
		T_neu = ph.refineAlg(max, twparameter);

		neu = T_neu.size();

		schrittweite = ph.SCHRITTWEITE;
		fehler = ph.FEHLER;*/
	}

	end = clock();

	cout << end - start << endl;

	delete viewer;

	return 0;
	

	int k = 0;

	while (max > tol) {

		if (viewer) viewer->CloseAll();


		AlpRiderPhase ph(neu);

		ph.newGrid(T_neu);

		ph.SCHRITTWEITE = schrittweite;
		ph.FEHLER = fehler;

		TWfolder folder(&twparameter, 0);
		folder.Add(&ph);
		folder.Init();


		//folder.Init(viewer);


		ph.FromMATLAB("alpTEST.m");

		folder.Loop(0);

		ph.ToMATLAB("alpTEST.m");

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
		AlpRiderPhase ph(neu);

		ph.newGrid(T_neu);

		ph.SCHRITTWEITE = schrittweite;
		ph.FEHLER = fehler;

		TWfolder folder(&twparameter, 0);
		folder.Add(&ph);
		folder.Init();

		folder.Init(viewer);

		ph.FromMATLAB("alpTEST.m");

		folder.Loop(0, 1);

		viewer->CloseAll();
	}

	end = clock();

	cout << end - start << "  -  Schritte:" << k + 1 << endl;

	delete viewer;

	return 0;
}
