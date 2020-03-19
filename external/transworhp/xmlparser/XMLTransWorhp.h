/*----------------------------------------------------------------
*
* Header-Datei fuer XMLTransWorhp.cpp
*
*----------------------------------------------------------------*/

#include <vector>
#include <sstream>
#include "TransWORHP.h"
#include "Variable.h"
#include "math.h"

using namespace std;

class XMLTransWorhp : public TransWorhp {

private:
	vector<Variable*> states; 	// Zustaende
	vector<Variable*> controls; // Steuerungen
	vector<Variable*> params; 	// Parameter
        vector<Variable*>::size_type stateCount, controlCount, paramCount, randCount, nebenCount;
	vector<Variable*> randc; 	// Randbedingungen
	vector<Variable*> nebenc; 	// Nebenbedingungen
	Node *objAnchor;			// Zielfunktion
	int **diffArray; 			// Ableitungs-Struktur

public:
	// Konstruktor
	XMLTransWorhp(int dis, vector<Variable*> _states, vector<Variable*> _controls, vector<Variable*> _params, Node* _objAnchor, vector<Variable*> _randc, vector<Variable*> _nebenc);

	// Hilfsfunktion stellt diffArray auf
	void makeStructure(Node *n, int num);

	// Zielfunktion
	double obj();
	double evaluateObj(Node *n);				// liest den Zielfunktionsbaum aus
	double evaluateDiffObj(Node *n, int diff_param); 	// leitet Zielfunktion ab
	bool obj_structure(DiffStructure &s);
	bool obj_diff(DiffStructure &s);

	// ODE-System
	void ode(double *dx, double t, const double *x, const double *u, const double *p);
	double evaluateOde(Node *n, const double *x, const double *u, const double *p); 			// liest ODE aus
	double evaluateDiffOde(Node *n, const double *x, const double *u, const double *p, int diff_param); 	// leitet ODE ab
	bool ode_structure(DiffStructure &s);
	bool ode_diff(DiffStructure &s, double t, const double *x, const double *u, const double *p);
	bool ode_diff_p(DiffStructure &s, double t, const double *x, const double *u, const double *p, int index);

	// Box Beschraenkungen	
	void x_boundary(double *x_low, double *x_upp);
	void u_boundary(double *u_low, double *u_upp);
	void p_boundary(double *p_low, double *p_upp);
	
	// Anfangs- und Endwerte
	void var_boundary(double *x_low, double *x_upp);

	// Randwerte
	void rand(double *r);
	void rand_boundary(double *r_low, double *r_upp);
	bool rand_structure(DiffStructure &s);
	bool rand_diff(DiffStructure &s);

	// Nebenbedingungen
	void neben(double *c, double t, const double *x, const double *u, const double *p);
	void neben_boundary(double *c_low, double *c_upp);
	bool neben_structure(DiffStructure &s);
	bool neben_diff(DiffStructure &s, double t, const double *x, const double *u, const double *p);
	bool neben_diff_p(DiffStructure &s, double t, const double *x, const
	double *u, const double *p, int index);

	// Startschaetzung
	void x_init(double *x, int i, int dis);
	void u_init(double *u, int i, int dis);
	void p_init(double *p);
};
