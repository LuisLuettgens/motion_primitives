#pragma once

#ifndef _WIN32
#include "TWGUIconfig.h"
#endif

#include "worhp/worhp.h"

#include "TWparameter.h"
#include "timing.h"
#include "diffstructure.h"

#ifdef TRANSWORHP_GRAPHICS
#include "Viewer.h"
#else

#include <string>

#ifdef WIN32
#ifdef _MSC_VER
#define DllExport __declspec( dllexport )
#pragma warning (disable : 4251)
#else
#define DllExport
#endif
#else
#define DllExport
#endif

#include "TransWORHP.h"

namespace tw {

class TransWorhpProblem;
class TransWorhp;

class DllExport Viewer {
public:
	Viewer(TWparameter* /*p*/) {}
	void init(TransWorhp */*p*/) {}
	void Matrix(const std::string &/*s*/, WorhpMatrix */*m*/) {}
	void Loop(int &) {}
	void update() {}
	int startTime, running;
	void autoScale() {}
	void closeAll() {}
};

int SDL_GetTicks();

struct FakeThread {
	int Active() {return 0;}
};
extern FakeThread *thethread0;

}

#endif

namespace tw {

/** Sammelordner fuer (verknuepfte) OCPs */
class DllExport TWfolder {

private:

	/** Flag ob Gitteranpassung aktiviert ist */
	bool meshRefFlag;

	/** Flag fuer Init */
	bool initFlag;

public:
	TWfolder(TWparameter *twparam, int n_con);
	virtual ~TWfolder();

	int Init(bool verbose=true, Params *params=nullptr);
	int Init(Viewer *view);

	int Reinit(bool verbose=0);

	void Add(TransWorhpProblem *t);

	/** Methode zur Gitteranpassung.
	* @author Matthias Rick
	* @return Fehlercode
	*/
	int meshRef();

	/** Methode zum testen der Loesung.
	* der Startwert wird hochintegriert und mit dem Endwert verglichen
	* @return hochintegrierter Wert und opt. Wert stimmen ueberein
	*/
	bool verify(bool verbose=false) const;

	/** Optimierungsaufruf */
	virtual int Loop(int wait=1, int terminate=0);
	virtual int WorhpLoop();

	void Show();

	void HM_structure(WorhpMatrix &DF, WorhpMatrix &DG, WorhpMatrix &HM);

	void HM_calculate();

	virtual bool step();

	/**
	* gibt an, welche Beschraenkung an der Stelle i ist.
	* @param i Kompomente des G Vektor
	* @return Typ/Index con
	*/
	std::string type_G(int i) const;

	std::vector<TransWorhpProblem*> phases;

	/** speichert: welche Phase beginnt in welcher Zeile */
	std::vector<int> phasenOffset;


	TWparameter *twparameter;

	bool verbose;


	bool Interrupt;

	OptVar worhp_o;
	Workspace worhp_w;
	Params worhp_p;
	Control worhp_c;

	Viewer *viewer;


	/** Start PHASE */

	/** Constrain-Matrix inititialisieren. */
	void init0(int nn, int mm);

	void GetSpace(int delta2); // Const. Matrix setzen



	// Spaltenweise DG-Struktur
	void DG_structure(int &ind, int colindex, WorhpMatrix *DG=0) const;

	/** Verbindung mit den Variablen von WORHP */
	void Connect();

	/** Aufruf der Randwerte (im Folder) */
	void Boundary();


	void Constraints();
	void Constraints2(double *GG);

	/** Zielfunktion.
	* bei mehreren Phasen: Summe der Zielfunktionen
	* alternativ: diese Methode ueberschreiben
	*/
	virtual double obj();

	/** Benutzerdefinierte Boundary
	* !nur linear moeglich!
	*/
	virtual void g_boundary(double *g_low, double *g_upp);

	/** Beschraenkungen.
	* zum Verbinden mehrerer Phasen
	*/
	virtual void con(double *GG);

	/** Struktur der Beschraenkungen.
	* optional
	*/
	virtual bool con_structure(DiffStructure &s);
	/** Ableitung der Beschraenkungen.
	* optional
	*/
	virtual bool con_diff(DiffStructure &s, int colindex);


	virtual void DG_calculate(int DG_start, int &ind, int colindex, WorhpMatrix &DG);

	/** Ableitungsstruktur */
	DiffStructure DS_con;

	double *G;
	double *G_low, *G_upp;
	double *Mu;

	int G_offset;
	int n_con;

	/** End PHASE */

	void WriteSensitivity(std::ostream &os, char var, char pert);
	void WriteSensitivityBinary(std::ostream &os, char var, char pert);
	void GetSensitivity(char var, char pert, int i, double *d);

	void PrintTimes(std::ostream &os);

	WORHP_Timing worhp_timing;
protected:

	std::vector<double> tmpG1;
	std::vector<double> tmpG2;

private:
	bool showOnly;

};

}
