/*----------------------------------------------------------------
 *
 * Einlesen von XML-Dateien
 *
 *----------------------------------------------------------------*/

#ifdef WIN32
#include "windows.h"
#endif

#include "XMLTransWorhp.h"

// Einlesefunktionen
Variable* setInfo(XMLNode *n, int isState);
void readFile(string name);	

// Variablen speichern alle Daten aus XMLFile
vector<Variable*> states;
vector<Variable*> controls;
vector<Variable*> params;
vector<Variable*> randc;
vector<Variable*> nebenc;
Node *objAnchor = new Node("");
double runtime = 1;

/////////////////////////////////////////////////////////////////////////////

// Funktion: parst XML Datei
void readFile(string name) {

	// Öffnen von XML-Dateien und Fehlerausgabe
	string datei = name;

	XMLParser p;
	
	XMLNode *node = p.Parse(datei);
	if (!node) {
		cout << "Error: file not found" << endl;
		exit (EXIT_FAILURE);	
	}

	// Einlesen der Daten
	XMLNode *mainnode = p.GetDocument();
	
	if (mainnode) {
		// Zielfunktion
		XMLNode *n = mainnode->GetFirstChild("objective");
		if (n) {
			objAnchor->setValue(n->GetText());
		}

		// Laufzeit
		n = mainnode->GetFirstChild("runtime");
		if (n) {
			runtime = atof(n->GetAttribute("end").c_str());
		}

		// Zustaende
		n = mainnode->GetFirstChild("state"); 
		while (n) {
			Variable *tmp = setInfo(n,1);
			states.push_back(tmp);
			n = mainnode->GetNextChild("state");
		}

		// Steuerungen
		n = mainnode->GetFirstChild("control");
		while (n) {
			Variable *tmp = setInfo(n,0);
			controls.push_back(tmp);			
			n = mainnode->GetNextChild("control");
		}
	
		// Parameter
		n = mainnode->GetFirstChild("param");
		while (n) {
			Variable *tmp = setInfo(n,0);
			params.push_back(tmp);			
			n = mainnode->GetNextChild("param");
		}

		// Randbedingungen
		n = mainnode->GetFirstChild("rand");
		while (n) {
			Variable *tmp = setInfo(n,2);
			randc.push_back(tmp);			
			n = mainnode->GetNextChild("rand");
		}

		// Nebenbedingungen
		n = mainnode->GetFirstChild("neben");
		while (n) {
			Variable *tmp = setInfo(n,2);
			nebenc.push_back(tmp);			
			n = mainnode->GetNextChild("neben");
		}
	}

	// Array mit den Variablen Namen erstellen, um Knoteninhalt zu interpretieren
	int stateCount = states.size();
	int controlCount = controls.size();
	int paramCount = params.size();

	string* stateNames = new string[stateCount];	
	string* controlNames = new string[controlCount];
	string* paramNames = new string[paramCount];

	vector<Variable*>::size_type i;

	for (i=0; i<stateCount; i++) 
		stateNames[i] = states[i]->getName();
	
	for (i=0; i<controlCount; i++) 
		controlNames[i] = controls[i]->getName();
	
	for (i=0; i<paramCount; i++) 
		paramNames[i] = params[i]->getName();
	
	// Baum aus Zielfunktion erstellen
	objAnchor->buildTree(stateNames, controlNames, paramNames, stateCount, controlCount, paramCount);

	// Baeume aus den DGL erstellen
	for (vector<Variable*>::size_type i=0; i < stateCount; i++) {
		states[i]->getTree()->buildTree(stateNames, controlNames, paramNames, stateCount, controlCount, paramCount);
	}

	// Baeume fuer die Randbedingungen erstellen
	for (vector<Variable*>::size_type i=0; i < randc.size(); i++) {
		randc[i]->getTree()->buildTree(stateNames, controlNames, paramNames, stateCount, controlCount, paramCount);
	}
	// Baeume fuer die Nebenbedingungen erstellen
	for (vector<Variable*>::size_type i=0; i < nebenc.size(); i++) {
		nebenc[i]->getTree()->buildTree(stateNames, controlNames, paramNames, stateCount, controlCount, paramCount);
	}
}

// Funktion: liest Werte aus XML Datei ein und fuellt damit ein Variable-Objekt
Variable* setInfo(XMLNode *n, int varNum) {
	// Erstellen eines neuen Variable-Objektes	
	Variable *var = new Variable(varNum);
	XMLNode *n1;
	
	if (varNum == 0 || varNum == 1) {
		var->setName(n->GetAttribute("name"));
	}

	n1 = n->GetFirstChild("boundaries");
	if (n1) {
		if (n1->GetAttribute("lower").compare("")) {
			var->setLowBound(atof(n1->GetAttribute("lower").c_str()));
		}
		if (n1->GetAttribute("upper").compare("")) {
			var->setUppBound(atof(n1->GetAttribute("upper").c_str()));
		}
	}

	// nur Zustaende und Rand- und Nebenbedingungen
	if (varNum == 1 || varNum == 2) {
		n1 = n->GetFirstChild("term");
		if (n1) {
			Node *anchor = new Node(n1->GetText());
			var->setTree(anchor);
		}

		// nur Zustaende
		if (varNum == 1) {
			n1 = n->GetFirstChild("start");
			if (n1) {
				var->setStart(atof(n1->GetText().c_str()));
			}

			n1 = n->GetFirstChild("final");
			if (n1) {
				var->setFinal(atof(n1->GetText().c_str()));
			}
		}
	}

	// Startschaetzung nur bei Zustaenden, Steuerungen und Parametern setzbar
	if (varNum == 1 || varNum == 0) {
		n1 = n->GetFirstChild("initialguess");
		if (n1) {
			var->setGuess(atof(n1->GetText().c_str()));
		}
	}

	return var;
}

/////////////////////////////////////////////////////////////////////////////

// Funktion: main-Funktion
int main(int argv, char* argc[]) {
	if (argv == 1) {
		cout << "Error: XML file not found" << endl;
		exit (EXIT_FAILURE);
	}

	readFile(argc[1]);

	TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv,argc);
	TWfolder folder(&twparameter,0);

	XMLTransWorhp *ph = new XMLTransWorhp(twparameter.NDIS,states,controls,params,objAnchor,randc,nebenc);

	// Default: LinearTimeAxis(0,1)
	if (runtime != 1) {
		ph->LinearTimeAxis(0,runtime);
	}

	folder.Add(ph);

	Viewer *viewer = 0;
	if (twparameter.PLOT) viewer = new Viewer(&twparameter);

	folder.Init();
	folder.Init(viewer);
	folder.Loop();

	delete viewer;

	cout << "\033[0m";
	return 0;
}