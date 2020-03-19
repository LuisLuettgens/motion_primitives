#pragma once

#ifndef _WIN32
#include "TWGUIconfig.h"
#endif

#include "xmlio.h"
#include "TWconsole.h"
#include "TransWORHP.h"
#include "LobattoPmTransWORHP.h"

#include <memory>

#ifdef TRANSWORHP_GRAPHICS
#include "../gui/sdlscreen.h"
#else
struct TWwindow {
	TWwindow() {}
};
#endif

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

/** Unterscheiden der Diskretisierung fuer TransWORHP */
class DllExport TWdiscretization {
public:

	// Standard = Trapez
	TWdiscretization(TWdiscretizationType t=TWdiscretizationType::Trapez, int s=1, int i=0) : type(t), stufen(s), innenpunkt(i) {}
	TWdiscretization(TWdiscretization& other) : type(other.type), stufen(other.stufen), innenpunkt(other.innenpunkt) {}

	/** Anzahl der Punkte der Intervalle 0 bis n_dis */
	int stuetzstellen(int n_dis) const {
		return (1+innenpunkt)*(n_dis-1)+1;
	}

	/** Anzahl der Punkte pro Intervall */
	int punkte() const {
		return 1+innenpunkt;
	}

	TWdiscretizationType type;

	int stufen;
	int innenpunkt;
};




class DllExport TWparameter {
public:
	TWparameter(std::string filename);
	~TWparameter();

	std::map<std::string,std::string>  Arguments(int argv, char* argc[]);

	void parseXML(const std::unique_ptr<XMLNode> &xmlmain);

	static std::unique_ptr<XMLNode> ReadParams(const std::string &filename);

	void printParams() const;

	/** Anzahl aller Parameter die eingelesen werden koennen */
	int countParams;
	/** Anzahl gesetzter/eingelesener Parameter */
	int setParams;

	/** Loesungsverfahren: volle Diskretisierung, Mehrfach-Schiessen, Pseudospektral */
	TransWORHP_type solver;

	/** Info ueber aktuelles Diskretisierungsschema */
	TWdiscretization twdiscretization;

	int butchertableau;
	double stepsize;
	double abserr, relerr;
	bool linInter; /**< lineare Interpolation der Steuerung: true=an, false=konst Steuer */
	bool parallel; /**< paralleles Intergrieren (explTW) mit openMP */

	/** Berechnungsart der Hessematrix */
	int hessianvalues;
	/** Struktur der Hessematrix */
	int hessianstructure;


	std::string paramfile_worhp;
	std::string paramfile_transworhp;
	int USERDF, USERDG, USERHM;

	bool showDF, showDG, showHM;
	bool showGrid;

	/** Struktur der Zielfunktion in Konsole ausgeben */
	bool showOBJstructure;
	/** Struktur der ODE in Konsole ausgeben */
	bool showODEstructure;
	/** Struktur der Randbedingungen in Konsole ausgeben */
	bool showRANDstructure;
	/** Struktur der Nebenbedingungen in Konsole ausgeben */
	bool showNEBENstructure;

	double eps;

	TWconsole twconsole;
	TWwindow twwindow;

	int NDIS;
	int PLOT;

	// Parameter Pseudospektral
	bool pm_smoothmode; // Polynom plotten
	int pm_displaypoints; // Ausfloesung fuer Plot
	PMnodeType pm_nodes_type; // welche Art von Stuetzstellen (Chebychev-Lobatto,Chebychev-Maxima, Legendre-Lobatto)

	// Parameter Mesh Refinement
	int meshref_mod; // Modus
	int meshref_err_mod; // Fehlerberuchnung Modus
	int meshref_M1; // Betts: M_1 Anz. neuer Punkte pro Intervall
	int meshref_R; //Betts: Ordnungsreduktion
	int meshref_M; //Betts: Anz. neuer Punkte Insgesamt
	int meshref_maxIter; // max. Anz. an Gitteranpassungsschritten
	double meshref_kappa; // Betts: kappa
	double meshref_tol; // Fehlerschranke

	int meshref_VER; // Zwischenschritte der Verfeinerung

	int meshref_PLOT_SW; //Plot Schrittweite
	int meshref_PLOT_ERR; //Plot Fehler
	int meshref_PLOT_MP; //Plot Position der Gitterpunkte
	int meshref_PLOT_LAMBDA; //Plot Beschraenkungen

	std::unique_ptr<XMLNode> xml;

	std::vector<int> multinode;
};

}
