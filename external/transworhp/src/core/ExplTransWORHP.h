#pragma once

#include "TransWORHP.h"

namespace tw {

/** struct zum Speichern von Informationen ueber die die Opt.Variablen */
struct typeX {
	int type; // 0: Zustand, 1: Steuerung, 2: Parameter
	int n; // Nummer
	int multinodeL; // Multiknoten links
	int dis; // diskreter Punkt
};

class DllExport ExplTransWorhp : public TransWorhp {

public:
	
	
	ExplTransWorhp(TransWorhpProblem *ph, TWparameter *twparam, std::vector<int> multinodes, std::vector<int> boxNB);
	
	/**
	*
	*/
	virtual ~ExplTransWorhp();
	
	
	
	int DoubleFrom(TransWorhp *ph) override;

	/** Hochintegrieren */
	int Integrate(int btableau) override;
	

private:
	int integrate_schrittweite(int index, int DGflag=0);
//#ifdef _OPENMP
	int integrate_schrittweite_parallel(int multi, int index, int DGflag=0);

	/** Blockweise tmp-Speicher fuer parallele Integration */
	std::vector<double> tmp_ode_parallel;
	/** Blockweise tmp-Speicher fuer parallele Integration */
	std::vector<double> tmp_ctrl_1_parallel;
	/** Blockweise tmp-Speicher fuer parallele Integration */
	std::vector<double> tmp_ctrl_2_parallel;
	
	/** Anzahl der Threads */
	int n_threads;
	
	/** Butcher-Objekte fuer parallele Integration */
	std::vector<Butcher> butcher_parallel;
//#endif
	std::vector<std::vector<double> > schrittweiten;

	
public:
	/** initialisiert das Butcher-Tableau.
	* @param index Tableau-Nummer
	* @param stepsize Schrittweite
	*/
	void butcherInit(int index, double stepsize=1e-3, bool verbose=false);
	
	/**
	* gibt an, welche Beschraenkung an der Stelle i ist (ODE,RAND,NB,etc).
	* @param Kompomente des G Vektor
	* @return Typ, zB ODE 0, RAND 1, etc
	*/
	std::string type_G(int i) const override;

	/** Zugriffsfunktionen */
	double x(int dis, int ode) const override;
	double u(int dis, int ctrl) const override;
	double p(int param) const override;

	/** Varianten fuer Zugriff */
	double x__(int dis, int ode) const override;
	double u__(int dis, int ctrl) const override;

	/** Indexbestimmung */
	int x_index(int dis, int ode) const override;
	int x_index__(int dis, int ode) const override;
	int u_index(int dis, int ctrl) const override;
	int u_index__(int dis, int ode) const override;
	int p_index(int param) const override;
	
	int x_indexode(int ode) const override;
	int u_indexode(int ctrl) const override;
	int p_indexode(int param) const override;
	
	/** laeufte ueber die Stuetzstellen, ausser die Multi-Knoten */
	int dis_index(int dis);
private:
	/** Knotenindezies ohne Multiknoten */
	std::vector<int> T_ohneMulti;
	/** speichert, ob Knoten ein Multiknoten ist.
	* nicht direkt drauf zugreifen, sondern ueber "isMultinote()"
	*/
	std::vector<bool> isMulti;
	
	/** speichert den Typ, die Nummer, dis. Punkt, linken Multiknoten der Opt.Varaiblen.
	* wird in init0() gefuellt
	*/
	std::vector<typeX> typeOfX;
public:
	/** @return 1, wenn Multiknoten, sonst 0 */
	bool isMultinode(int i);
	
	/** fuer ROW Verfahren muss diese Methode implementiert werden */
	virtual bool ode_diff2(DiffStructure &ds, double t, const double *x, const double *u, const double *p);

	void GetState(double *x, double t) override;
	void GetControl(double *uu, double t) override;
	
#ifdef TRANSWORHP_GRAPHICS
	std::string setViewerTitle() override;
	void updateViewer(DataStorage* ds, std::vector<double> &tmptime) override;
	void setTemptimeForViewer(std::vector<double> &tmptime) override;
#endif
	
	void GetBoundaryIndices(std::vector<int> &indices, int d) override;
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

public:
	double Objective(double ScaleObj) override;

	int DF_structure(WorhpMatrix *DF=nullptr, int offset=0) override;
	void DF_calculate(WorhpMatrix &DF, double ScaleObj) override;

	int DG_structure(const TWfolder *f, WorhpMatrix *DG=nullptr, int offset=0) override;
	void DG_calculate(TWfolder *f, WorhpMatrix &DG) override;
	
	int HM_structure_ohne_Diag(int, WorhpMatrix &DF, WorhpMatrix &DG, WorhpMatrix *HM=nullptr, int offset=0) override;
	void HM_calculate1(WorhpMatrix &DF, WorhpMatrix &DG, WorhpMatrix &HM, double ScaleObj, double *Mu) override;
	
	void Lagrange() override;

	void init0() override;

	/** Mehrzielknoten. */
	std::vector<int> multinodes;
	int n_multinodes;

	/* Zugriffsmethoden */
	/** Zugriff auf Zustaende.
	* speichert auch temp. Zustaende, die keine Optimierungsvariablen sind
	*/
	std::vector<double*> pointer_x;
	std::vector<double*> pointer_u;
	double *pointer_p;

protected:

	/** Stetigkeitsbedingung fuer Mehrfachschiessen. */
	void Continuous(double* ddx, int dis);

	std::vector<double> block_x;
	std::vector<double> block_u;

	/** Hochintegrieren der Zustaende */
	void IntegrateStates();
	/** Hochintegrieren der Zustaende (parallel) */
	void IntegrateStates_parallel();
	/** Hochintegrieren der Zustaende in Intervall [start,end] */
	void IntegrateStates2(int start, int end, int DGflag=0);
	/** Hochintegrieren der Zustaende als Startschaetzung.
	* bis auf das setzen der Werte an den Multiknoten identisch zu "IntegrateStates"
	* @return Anzahl der Integrationsschritte
	*/
	int IntegrateInitial();

	void fromMATLAB_impl(std::ifstream &stream) override;

private:
	/** beschraenkte Zustaende */
	const std::vector<int> boxNB;

public:
	/** Anzahl der Zustandesbeschraenkungen */
	int n_boxNeben;

private:
	/* Methoden zum Erweitern der Box-Schranken auf Nebenbedingungen
	fuer Zustande die keine Optimierungsvariable sind */
	void boxNeben(double *c, const double *x);
	void boxNeben_boundary(double *c_low, double *c_upp);
};

}
