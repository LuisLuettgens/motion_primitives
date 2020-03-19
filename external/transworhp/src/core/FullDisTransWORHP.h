#pragma once

#include "TransWORHP.h"

namespace tw {

class XMLNode;
class TWbaseSpline;

/**
* struct fuer (internen) Export von TW Daten
* @author Matthias Rick
*/
struct exportTW {
	std::vector<double> T; // Zeit
	std::vector<double> p; // Parameter
	std::vector<std::vector<double> > x; // Zustaende
	std::vector<std::vector<double> > u; // Steuerung
	std::vector<std::vector<double> > La; // Lambda
	std::vector<std::vector<double> > Mu; // Mu
	TWdiscretization *twdiscretisation; // Diskretisierung
	bool ok; // alles richtig gesetzt?
};


/** TransWorhp */
class DllExport FullDisTransWorhp : public TransWorhp {

public:

	FullDisTransWorhp(TransWorhpProblem *ph, TWparameter *twparam);
	
	/**
	* Destruktor
	*/
	virtual ~FullDisTransWorhp();
	
	/**
	* exportiert Lambda in m-Datei.
	* @param filename Dateiname
	*/
	void ToMATLAB_LambdaMu(const std::string& filename); //Matthias Rick
	/**
	* importiert Lambda aus m-Datei.
	* @param filename Dateiname
	*/
	void FromMATLAB_LambdaMu(const std::string& filename); //Matthias Rick

	int Integrate(int btableau) override;
	int integrate(int index); //TODO besserer Name? verwirrrt doch arg.
	int integrateRKF(int index, int &startflag);


	virtual void Debug_G();
	virtual void Debug_G2();

	/**
	* gibt an, welche Beschraenkung an der Stelle i ist (ODE,RAND,NB,etc).
	* @param Kompomente des G Vektor
	* @return Typ, zB ODE 0, RAND 1, etc
	*/
	std::string type_G(int i) const override;

	friend std::ostream& operator<<(std::ostream &os, const TransWorhp &p);
	

	/* --------------- Zugriffsfunktionen --------------- */
	/** gibt den Zustand an Gitterpunkt zurueck.
	* @param dis dis. Punkt
	* @param ode Zustand
	* @return Zustand 'ode' an Stelle 'dis'
	*/
	double x(int dis, int ode) const override;
	/** gibt die Steuerung an Gitterpunkt zurueck.
	* @param dis dis. Punkt
	* @param ctrl Steuerung
	* @return Steuerung 'ctrl' an Stelle 'dis'
	*/
	double u(int dis, int ctrl) const override;
	/** gibt freien Parameter zurueck.
	* @param param freier Parameter
	* @return freier Parameter
	*/
	double p(int param) const override;

	/* --------------- Varianten fuer Zugriff --------------- */
	/** gibt den Zustand an Gitterpunkt zurueck.
	* gibt den Zustand an Gitterpunkt zurueck, wobei auch Zwischenpunkte
	* beachtet werden (zB bei HermiteSimpson)
	* @param dis dis. Punkt
	* @param ode Zustand
	* @return Zustand 'ode' an Stelle 'dis'
	*/
	double x__(int dis, int ode) const override;
	/** gibt die Steuerung an Gitterpunkt zurueck.
	* gibt die Steuerung an Gitterpunkt zurueck, wobei auch Zwischenpunkte
	* beachtet werden (zB bei HermiteSimpson)
	* @param dis dis. Punkt
	* @param ctrl Steuerung
	* @return Steuerung 'ctrl' an Stelle 'dis'
	*/
	double u__(int dis, int ctrl) const override;

	/* --------------- Indexbestimmung --------------- */
	/** gibt den Index im Optimierungsvektor zurueck.
	* @param dis dis. Punkt
	* @param ode Zustand
	* @return Index
	*/
	int x_index(int dis, int ode) const override;
	/** gibt den Index im Optimierungsvektor zurueck.
	* es werden auch Zwischenpunkte beruecksichtigt
	* @param dis dis. Punkt
	* @param ode Zustand
	* @return Index
	*/
	int x_index__(int dis, int ode) const override;
	/** gibt den Index im Optimierungsvektor zurueck.
	* @param dis dis. Punkt
	* @param ctrl Steuerung
	* @return Index
	*/
	int u_index(int dis, int ctrl) const override;
	/** gibt den Index im Optimierungsvektor zurueck.
	* es werden auch Zwischenpunkte beruecksichtigt
	* @param dis dis. Punkt
	* @param ctrl Steuerung
	* @return Index
	*/
	int u_index__(int dis, int ctrl) const override;
	/** gibt den Index im Optimierungsvektor zurueck.
	* @param param freier Parameter
	* @return Index
	*/
	int p_index(int param) const override;

	/** gibt den Index des Zustandes zureuck.
	* @param ode Zustand
	* @return Index
	*/
	int x_indexode(int ode) const override;
	/** gibt den Index der Steuerung zureuck.
	* @param ctrl Steuerung
	* @return Index
	*/
	int u_indexode(int ctrl) const override;
	/** gibt den Index der freien Parameter zureuck.
	* @param param freier Paramter
	* @return Index
	*/
	int p_indexode(int param) const override;
	
	/** Grafische Ausgabe */
#ifdef TRANSWORHP_GRAPHICS
	
	/* set title of Viewer */
	std::string setViewerTitle() override;
	
	/* set state, control, etc for Viewer */
	void updateViewer(DataStorage *ds, std::vector<double> &tmptime) override;
	
	/* allocates memory and sets temptime for Viewer */
	void setTemptimeForViewer(std::vector<double> &tmptime) override;
#endif
	
	void GetState(double *x, double t) override;
	void GetControl(double *uu, double t) override;
	
	
protected:
	double &setx(int dis, int ode) override;
public:
	/** Speicher und Struktur */
	void GetSpace(int delta1, int delta2) override;
	/** Verbindung zu WORHP */
	void Connect(const OptVar &o, const Params &p) override;
	/** Beschraenkungsschranken */
	void Boundary() override;
	/** Beschraenkungen mit anderem G Vektor */
	void Constraints2(double *GG, int DGflag=0) override;
	
	double Objective(double ScaleObj) override;


protected:
	virtual void DG_diff_ode(double t, int k, int active_index);
	virtual void DG_diff_ode_p(double t, int l, int dis, int active_index);

	virtual void DG_diff_rand(int k);

	virtual void DG_diff_neben(double t, int dis);
	virtual void DG_diff_neben_p(double t, int l, int dis);

public:

	int DF_structure(WorhpMatrix *DF=nullptr, int offset=0) override;
	void DF_calculate(WorhpMatrix &DF, double ScaleObj) override;

	int DG_structure(const TWfolder *f, WorhpMatrix *DG=nullptr, int offset=0) override;
	void DG_calculate(TWfolder *f, WorhpMatrix &DG) override;
	
	int HM_structure_ohne_Diag(int hessianstructure, WorhpMatrix &DF, WorhpMatrix &DG, WorhpMatrix *HM=nullptr, int offset=0) override;
	void HM_calculate1(WorhpMatrix &DF, WorhpMatrix &DG, WorhpMatrix &HM, double ScaleObj, double *Mu) override;
	void HM_calculate2(WorhpMatrix &DF, WorhpMatrix &DG, WorhpMatrix &HM, double ScaleObj, double *Mu) override;
	void HM_calculate3(WorhpMatrix &DF, WorhpMatrix &DG, WorhpMatrix &HM, double ScaleObj, double *Mu) override;
	void HM_calculate4(WorhpMatrix &DF, WorhpMatrix &DG, WorhpMatrix &HM, double ScaleObj, double *Mu) override;
	void HM_calculate5(WorhpMatrix &DF, WorhpMatrix &DG, WorhpMatrix &HM, double ScaleObj, double *Mu) override;

	void Lagrange() override;

	void init0() override;

	void GetBoundaryIndices(std::vector<int> &indices, int d) override;

protected:
	
	WorhpMatrix DFtemp1, DGtemp1, DFtemp2, DGtemp2;
	bool tempmatrixinitialised;

	void tempmatrixinit();

	virtual void HermiteSimpson(double* G1, double *G2, int dis);
	virtual void Lobatto(double* G1, double *G2, double *G3, int dis);
	virtual void Trapez(double* ddx, int dis);
	virtual void Euler(double* ddx, int dis);
	virtual void RechteSeite(double *G1, double t, int dis);

	void fromMATLAB_impl(std::ifstream &stream) override;


/*Methoden fuer Dis.fehler - Matthias Rick*/
private:
	/**
	* Wertet die Differenz von Dynamik und Spline an Punkt t aus (vgl. Betts)
	* -> Hilfsfunktion fuer disFehlerIntegral
	* @author Matthias Rick
	* @param t Zeitpunkt
	* @param komp Komponente
	* @param zustand Zustand
	* @param steuerung Steuerung
	* @return |Dynamik - Spline|
	*/
	double auxFehlerInt(double t, int komp, TWbaseSpline **zustand, TWbaseSpline **steuerung);
	/**
	* Berechnet Integral ( |spline(x)' - f(spline(x),spline(u),t)| )
	* -> Hilfsfunktion fuer diskretisierungsfehlerBetts
	* @author Matthias Rick
	* @param fehler Rueckgabewert: Fehler pro Intervall
	* @param a Intervall Anfang
	* @param b Intervall Ende
	* @param zustand Zustand
	* @param steuerung Steuerung
	*/
	void disFehlerIntegral(double *fehler,const double a, const double b, TWbaseSpline **zustand, TWbaseSpline **steuerung);
	/**
	* berechnet den Diskretisierungsfehler nach Betts
	* @author Matthias Rick
	*/
	void diskretisierungsfehlerBetts();
	/**
	* berechnet den Diskretisierungsfehler durch Vergleich mit Verfahren hoeherer Ordnung
	* Euler mit Trapez
	* Trapez mit Hermite-Simpson
	* @author Matthias Rick
	*/
	void diskretisierungsfehler2();
public:
	/** Diskretisierungsfehler */
	std::vector<std::vector<double> > FEHLER;
	/** Diskretisierung */
	std::vector<std::vector<double> > SCHRITTWEITE;

	/**
	* Berechnung des Diskretisierungsfehler.
	* berechnet den Diskretisierungsfehler und schreibt ihn in FEHLER
	* @author Matthias Rick
	*/
	void diskretisierungsfehler();

	/**
	* Bestimmt das neue Gitter
	* @author Matthias Rick
	* @param max aktuelles Maximum (groesster Fehler)
	* @return Vektor mit neuer Diskretisierung
	*/
	std::vector<double> refineAlg(double &max);
private:
	/**
	* Berechnet das neue Gitter nach Betts
	* @author Matthias Rick
	* @param I Vektor in den die Anzahl an neuen Punkten pro Intervall geschrieben wird
	* @param max aktuelles Maximum
	*/
	void betts(std::vector<int> &I, double &max);
	/**
	* Berechnet Auf-/Absprungpunkte nach Bueskens
	* @author Matthias Rick
	* @return Vektor mit neuen Punkten
	*/
	std::vector<double> bueskens();

public:
	
	/**
	* exportieren von TW Daten
	* Zeit, Param, Zustand, Steuer, Lambda, Mu
	* @author Matthias Rick
	* @return exportTW-Struct
	*/
	exportTW outTW();

	/**
	* importieren von TW Daten
	* Zeit, Param, Zustand, Steuer, Lambda, Mu
	* @author Matthias Rick
	* @param eTW TW-Struct
	*/
	void inTW(exportTW eTW);
};

}
