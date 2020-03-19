/*----------------------------------------------------------------
 *
 * Example: Splineproblem
 *
 *----------------------------------------------------------------*/

#ifdef WIN32
#include "windows.h"
#endif

#include "TransWORHP.h"

#include "TWspline.h"

// Zeitmessung
#include <time.h>

using namespace std;

class ChemiePhase : public TransWorhp {
public:

	ChemiePhase(int dis) : TransWorhp("Chemiereaktor", dis, 2, 1, 0, 0, 0) {

		/** �ffnen von XML-Dateien (und Fehlerausgabe, wenn's schief l�uft) */
		string a = "transworhp.xml";
		XMLParser p;
		TextOutputType_e tot = ASCII;

		XMLNode *node = p.Parse(a); 	// Parses the file and returns root node.

		// std::cout << *node;			// Wrong! Don't write output directly to std::cout.

		DebugStream d(tot); 			// Open DebugStream.
		d << beginfile;
		if (node) {
			//d << *node; 			    // Write whole file.
			//d << *p.GetDocument(); 	// Write only main part.

		}
		else {
			p.GetError(d);			// Show all errors.
			//d << *node;			// No! We don't have this node,
			//d << *p.GetDocument(); 	// and no document at all.
		}
		d << endfile;

		std::cout << d.GetString();		// Finally write DebugStream anywhere you want.


		/** Ab hier: anschauen der Daten... */
		XMLNode *mainnode = p.GetDocument();
		if (mainnode) {
			XMLNode *n = mainnode->GetFirstChild("WORHP"); // Knoten WORHP holen
			if (n) {
				cout << n->GetAttribute("param") << endl; // sollte worhp.xml ausgeben
			}

			n = mainnode->GetFirstChild("DISCRETIZATION");
			if (n) {
				cout << n->GetText() << endl; // sollte 1 ausgeben
			}

			n = mainnode->GetFirstChild("PLOT");
			if (n) {
				XMLNode *n2 = n->GetFirstChild("SPARSITY"); // eine Ebene tiefer...
				while (n2) {                                // alle Knoten durchlaufen
					cout << n2->GetText() << endl; // sollte DF,DG,HM ausgeben
					n2 = n->GetNextChild("SPARSITY");
				}
			}




		}


	}

	double obj() {

		return -x(n_dis - 1, 1);
	}


	bool obj_structure(DiffStructure &s) {

		s(0, x_index(n_dis - 1, 1));
		return true;

	}

	void ode(double *dx, double t, const double *x, const double *u, const double *p) {
		dx[0] = -u[0] * x[0] + u[0] * u[0] * x[1];
		dx[1] = u[0] * x[0] - 3 * u[0] * u[0] * x[1];
	}


	// Optional
	bool ode_structure(DiffStructure &s) {
		//return false;
		s(0, u_index(0, 0));
		s(0, x_index(0, 0));
		s(0, x_index(0, 1));


		s(1, u_index(0, 0));
		s(1, x_index(0, 0));
		s(1, x_index(0, 1));

		return true;
	}

	void u_boundary(double *u_low, double *u_upp) {

		u_low[0] = 0;
		u_upp[0] = +1;

	}


	void var_boundary(double *x_low, double *x_upp) {

		x_upp[x_index(0, 0)] = x_low[x_index(0, 0)] = 1;
		x_upp[x_index(0, 1)] = x_low[x_index(0, 1)] = 0;

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
		//viewer->restrictionPlot(T,Lambda,n_dis,n_ctrl,n_ode,twdiscretization);
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

	int alt;
	int neu;
	

	const double tol = twparameter.meshref_tol;
	double max = 1.0;

	Viewer *viewer = 0;

	if (twparameter.PLOT) viewer = new Viewer(&twparameter);


	// Startschaetzung in Datei erstellen
	{
		ChemiePhase ph(twparameter.NDIS);
		
		TWfolder folder(&twparameter, 0);
		folder.Add(&ph);
		folder.Init();
		folder.Init(viewer); //<------------------------- grafik ausgabe

		folder.Loop(0);
		

		ph.ToMATLAB("chemieTEST.m");
		ph.ToMATLAB_LambdaMu("chemieLM.m");
		
		T_alt = ph.getGrid();
		T_neu = ph.refineAlg(max);

		neu = T_neu.size();

		schrittweite = ph.SCHRITTWEITE;
		fehler = ph.FEHLER;
	}
	
	return 0; 

	int k = 0;

	while (max > tol) {

		if (viewer) viewer->CloseAll();


		ChemiePhase ph(neu);

		ph.newGrid(T_neu);

		ph.SCHRITTWEITE = schrittweite;
		ph.FEHLER = fehler;

		TWfolder folder(&twparameter, 0);
		folder.Add(&ph);
		folder.Init();

		
		//folder.Init(viewer);
		
		
		ph.FromMATLAB("chemieTEST.m");
		ph.FromMATLAB_LambdaMu("chemieLM.m");

		folder.Loop(0);

		ph.ToMATLAB("chemieTEST.m");
		ph.ToMATLAB_LambdaMu("chemieLM.m");

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
		ChemiePhase ph(neu);

		ph.newGrid(T_neu);

		ph.SCHRITTWEITE = schrittweite;
		ph.FEHLER = fehler;

		TWfolder folder(&twparameter, 0);
		folder.Add(&ph);
		folder.Init();

		
		folder.Init(viewer);
		

		ph.FromMATLAB("chemieTEST.m");
		ph.FromMATLAB_LambdaMu("chemieLM.m");
		
		folder.Loop(0,1);

	}

	end = clock();

	cout << end - start << "  -  Schritte:" << k + 1 << endl;


	delete viewer;


	/*
	TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv,argc);
	TWfolder folder(&twparameter,0);

	ChemiePhase ph(twparameter.NDIS);
	folder.Add(&ph);

	Viewer *viewer = 0;
	if (twparameter.PLOT) viewer = new Viewer(&twparameter);

	folder.Init();
	folder.Init(viewer);
	folder.Loop();


	delete viewer;
	*/
	return 0;
}

