#pragma once

#ifndef _WIN32
#include "TWGUIconfig.h"
#endif

#include "TWcount.h"
#include "diffstructure.h"
#include "butcher.h"

#include "worhp/worhp.h"

#ifndef DllExport
#ifdef _WIN32
#define DllExport __declspec(dllexport)
#else
#define DllExport
#endif
#endif

#ifdef _MSC_VER
#pragma warning(disable : 4251)
#endif

namespace tw {

class TransWorhpProblem;
class TWfolder;
class DataStorage;
class Viewer;
class TWdiscretization;
class TWparameter;

#include <vector>

#ifndef M_PI
#define M_PI 3.14152926535897932846
#endif

enum class TransWORHP_type {unknown=-1, fullDiscretization=0, multipleShooting=2, pseudospectral=4, pseudospectral_gauss};

void DllExport worhpoutput(int i, const char *message);

class DllExport TransWorhp {

public:

	/**
	 * OCP Loeser
	 *
	 */
	TransWorhp(TransWorhpProblem *ph, TWparameter *twparam);

	/**
	* Destruktor
	*/
	virtual ~TransWorhp();


	TransWorhpProblem *phase;

	static void MatrixToMATLAB(const WorhpMatrix &m, const std::string& filename);
	static void DoubleToMATLAB(double *data, int n, const std::string& filename);


	virtual int DoubleFrom(TransWorhp *ph);

	/**
	* exportiert Zustand und Steuerung in m-Datei.
	* @param filename Dateiname
	*/
	virtual void ToMATLAB(const std::string& filename) const;
	/**
	* importiert Zustand und Steuerung aus m-Datei.
	* ruft auf dem solver fromMATLAB_impl() auf
	* @param filename Dateiname
	*/
	void FromMATLAB(const std::string& filename);

	virtual int Integrate(int btableau) = 0;


	/**
	* gibt an, welche Beschraenkung an der Stelle i ist (ODE,RAND,NB,etc).
	* @param Kompomente des G Vektor
	* @return Typ, zB ODE 0, RAND 1, etc
	*/
	virtual std::string type_G(int i) const = 0;

	void Structure_Sizes(TWfolder *f, int hessianvalues, int &DF_nnz, int &DG_nnz, int &HM_nnz);
	virtual void LinearTimeAxis(double start, double end);


	/* --------------- Zugriffsfunktionen --------------- */
	/** gibt den Zustand an Gitterpunkt zurueck.
	* @param dis dis. Punkt
	* @param ode Zustand
	* @return Zustand 'ode' an Stelle 'dis'
	*/
	virtual double x(int dis, int ode) const = 0;
	/** gibt die Steuerung an Gitterpunkt zurueck.
	* @param dis dis. Punkt
	* @param ctrl Steuerung
	* @return Steuerung 'ctrl' an Stelle 'dis'
	*/
	virtual double u(int dis, int ctrl) const = 0;
	/** gibt freien Parameter zurueck.
	* @param param freier Parameter
	* @return freier Parameter
	*/
	virtual double p(int param) const = 0;

	/* --------------- Varianten fuer Zugriff --------------- */
	/** gibt den Zustand an Gitterpunkt zurueck.
	* gibt den Zustand an Gitterpunkt zurueck, wobei auch Zwischenpunkte
	* beachtet werden (zB bei HermiteSimpson)
	* @param dis dis. Punkt
	* @param ode Zustand
	* @return Zustand 'ode' an Stelle 'dis'
	*/
	virtual double x__(int dis, int ode) const = 0;
	/** gibt die Steuerung an Gitterpunkt zurueck.
	* gibt die Steuerung an Gitterpunkt zurueck, wobei auch Zwischenpunkte
	* beachtet werden (zB bei HermiteSimpson)
	* @param dis dis. Punkt
	* @param ctrl Steuerung
	* @return Steuerung 'ctrl' an Stelle 'dis'
	*/
	virtual double u__(int dis, int ctrl) const = 0;

	/* --------------- Indexbestimmung --------------- */
	/** gibt den Index im Optimierungsvektor zurueck.
	* @param dis dis. Punkt
	* @param ode Zustand
	* @return Index
	*/
	virtual int x_index(int dis, int ode) const = 0;
	/** gibt den Index im Optimierungsvektor zurueck.
	* es werden auch Zwischenpunkte beruecksichtigt
	* @param dis dis. Punkt
	* @param ode Zustand
	* @return Index
	*/
	virtual int x_index__(int dis, int ode) const = 0;
	/** gibt den Index im Optimierungsvektor zurueck.
	* @param dis dis. Punkt
	* @param ctrl Steuerung
	* @return Index
	*/
	virtual int u_index(int dis, int ctrl) const = 0;
	/** gibt den Index im Optimierungsvektor zurueck.
	* es werden auch Zwischenpunkte beruecksichtigt
	* @param dis dis. Punkt
	* @param ctrl Steuerung
	* @return Index
	*/
	virtual int u_index__(int dis, int ctrl) const = 0;
	/** gibt den Index im Optimierungsvektor zurueck.
	* @param param freier Parameter
	* @return Index
	*/
	virtual int p_index(int param) const = 0;

	/** gibt den Index des Zustandes zureuck.
	* @param ode Zustand
	* @return Index
	*/
	virtual int x_indexode(int ode) const = 0;
	/** gibt den Index der Steuerung zureuck.
	* @param ctrl Steuerung
	* @return Index
	*/
	virtual int u_indexode(int ctrl) const = 0;
	/** gibt den Index der freien Parameter zureuck.
	* @param param freier Paramter
	* @return Index
	*/
	virtual int p_indexode(int param) const = 0;

	/** Grafische Ausgabe */
#ifdef TRANSWORHP_GRAPHICS

	/* set title of Viewer */
	virtual std::string setViewerTitle() = 0;

	/* set state, control, etc for Viewer */
	virtual void updateViewer(DataStorage *ds, std::vector<double> &tmptime) = 0;

	/* allocates memory and sets temptime for Viewer */
	virtual void setTemptimeForViewer(std::vector<double> &tmptime) = 0;
#endif

	void TimeAxis(double exponent);

	virtual void GetState(double *x, double t) = 0;
	virtual void GetControl(double *uu, double t) = 0;

	virtual void PrintOCPstates(std::ostream *os=nullptr) const;
	void PrintMultipliers(std::ostream *os=nullptr) const;

protected:
	virtual double &setx(int dis, int ode) = 0;
public:
	/** Speicher und Struktur */
	virtual void GetSpace(int delta1, int delta2) = 0;
	/** Verbindung zu WORHP */
	virtual void Connect(const OptVar &o, const Params &p) = 0;
	/** Beschraenkungsschranken */
	virtual void Boundary() = 0;

protected:
	/** Prueft, ob obere und untere Schranke fuer OptVar gleich sind
	 * und passt dann den Startwert von X an, damit dieser zulaessig ist */
	void boxConToInitGuess();

	/** Gibt die Struktur von
	 * Obj,ODE,Rand,Neben
	 * aus */
	void printStructure() const;
public:
	/** Beschraenkungen */
	void Constraints();
	/** Beschraenkungen mit anderem G Vektor */
	virtual void Constraints2(double *GG, int DGflag=0) = 0;

	virtual double Objective(double ScaleObj) = 0;

	virtual int DF_structure(WorhpMatrix *DF=nullptr, int offset=0) = 0;
	virtual void DF_calculate(WorhpMatrix &DF, double ScaleObj) = 0;

	virtual int DG_structure(const TWfolder *f, WorhpMatrix *DG=nullptr, int offset=0) = 0;
	virtual void DG_calculate(TWfolder *f, WorhpMatrix &DG) = 0;

	int HM_structure(int hessianstructure, WorhpMatrix &DF, WorhpMatrix &DG, WorhpMatrix *HM=nullptr, int offset=0);
	virtual int HM_structure_ohne_Diag(int hessianstructure, WorhpMatrix &DF, WorhpMatrix &DG, WorhpMatrix *HM=nullptr, int offset=0) = 0;
	void HM_calculate(int hessianvalues, WorhpMatrix &DF, WorhpMatrix &DG, WorhpMatrix &HM, double ScaleObj, double *Mu);
	virtual void HM_calculate1(WorhpMatrix &DF, WorhpMatrix &DG, WorhpMatrix &HM, double ScaleObj, double *Mu) = 0;
	virtual void HM_calculate2(WorhpMatrix &DF, WorhpMatrix &DG, WorhpMatrix &HM, double ScaleObj, double *Mu);
	virtual void HM_calculate3(WorhpMatrix &DF, WorhpMatrix &DG, WorhpMatrix &HM, double ScaleObj, double *Mu);
	virtual void HM_calculate4(WorhpMatrix &DF, WorhpMatrix &DG, WorhpMatrix &HM, double ScaleObj, double *Mu);
	virtual void HM_calculate5(WorhpMatrix &DF, WorhpMatrix &DG, WorhpMatrix &HM, double ScaleObj, double *Mu);

	virtual void Lagrange() = 0;

	virtual void init0() = 0;

	virtual void GetBoundaryIndices(std::vector<int> &indices, int d) = 0;

protected:

	int getncon();
	int getnvar();
	
	virtual void fromMATLAB_impl(std::ifstream &stream) = 0;

public:
	/**
	* gibt das aktuelle Gitter zurueck
	* @author Matthias Rick
	* @return Gitterpunkte als Vektor
	*/
	std::vector<double> getGrid() const;
	/**
	* setzt das Zeitgitter
	* @param grid neues Gitter
	* @author Matthias Rick
	*/
	void newGrid(std::vector<double> grid);

	int &n_dis;
	const int &n_ode; /**< Anzahl Zustaende */
	const int &n_ctrl; /**< Anzahl Steuerungen */
	const int &n_param; /**< Anzahl freier Parameter */
	const int &n_rand; /**< Anzahl Randbedingungen */
	const int &n_neben; /**< Anzahl Nebenbedingungen */
	const int &n_integral; /**< Anzahl Integralterm */

	int Delta1; /**< Offset fuer WORHP X */
	int Delta2; /**< Offset fuer WORHP G */
	int n_var; /**< Anzahl Optimierungsvariablen */
	int n_con; /**< Anzahl Beschraenkungen */
	int n_zen;

	std::vector<double> T; /**< Zeit-Diskretisierung */

	double *X; /**< Optimierungsvektor */
	double *G; /**< Beschraenkungen */
	double *ZEN; // Zen

	double *X_low, *X_upp;
	double *G_low, *G_upp;

	double *Lambda, *Mu;

	Viewer *viewer; /**< Ausgabefenster */

	/** TransWorhp Parameter */
	TWparameter *twparameter;
	/** Info ueber aktuelles Diskretisierungsschema */
	TWdiscretization *twdiscretization;

	Butcher butcher; /**< Integrationsverfahren */

	/** Loesungsverfahren.
	* 0 = direktes (normales) Verfahren
	* 2 = Mehrzielmethode
	* 4 = Pseudospektral-Methode
	*/
	TransWORHP_type transworhp_type;

	TWcount_calls twcount_calls;

	int DF_start, DF_nnz;
	int DG_start, DG_nnz;
	int HM_start, HM_start_diag, HM_nnz_ohne_diag;

	double Infty; /**< WORHP-Unendlich */ //public: damit es von aussen ueberpruefbar ist



protected:
	int SHOWDF, SHOWDG, SHOWHM;

	const double eps;

	int USERDF, USERDG, USERHM;

	// Hilfsvariablen fuer Auswertung ODE, NB, etc.
	std::vector<double> tmp_ode_1, tmp_ode_2;
	std::vector<double> tmp_ctrl_1, tmp_ctrl_2;
	std::vector<double> tmp_rand_1, tmp_rand_2;
	std::vector<double> tmp_neben_1, tmp_neben_2;
	std::vector<double> tmp_integral_1, tmp_integral_2, tmp_integral_12;

	std::vector<double> tmp_gg_1, tmp_gg_2, tmp_gg_3, tmp_gg_4, tmp_gg_5;

	DiffStructure DS_obj, DS_ode, DS_neben, DS_rand, DS_integral;

public:
	std::vector<double> lagrange_integral;
	std::vector<double> lagrange_weight;
};

}
