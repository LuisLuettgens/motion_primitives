#pragma once

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

#ifndef _WIN32
#include "TWGUIconfig.h"
#endif

#ifdef TRANSWORHP_GRAPHICS
namespace tw {
class Viewer;
}
#endif

#include <memory>
#include <vector>

namespace tw {

class TWparameter;
class TransWorhp;
class DiffStructure;

struct DllExport TWdimension {

	std::string ID;
	int n_dis;
	int n_ode;
	int n_ctrl;
	int n_param;
	int n_rand;
	int n_neben;
	int n_integral;
	int n_zen;
	//only for explTW
	std::vector<int> multinode;
	std::vector<int> BOXneben;

	TWdimension();

	TWdimension(std::string name, int dis, int ode, int ctrl, int param, int rand, int neben,
	            int integral = 0, int zen = 0);

	/** setzt Mehrzielknoten
	 * @param multi Anzahl der Mehrzielknoten
	 */
	void setMultinodes(int multi);
};


class DllExport TransWorhpProblem {

public:

	TransWorhpProblem(const TWdimension &TWdata);

	/** default Destruktor */
	virtual ~TransWorhpProblem() = default;

	std::unique_ptr<TransWorhp> solver;

	void setSolver(TWparameter *twparam);


	/* --------------- Implementierung des OCP --------------- */

	/** Zielfunktion.
	* @return Zielfunktionswert
	*/
	virtual double obj()=0;
	/** Struktur der Zielfunktion.
	* optional
	* @param s Abhaengigkeit
	* @return true, wenn Struktur benutzt werden soll
	*/
	virtual bool obj_structure(DiffStructure &s);
	/** Ableitung der Zielfunktion.
	* optional
	* @param s Abhaengigkeit
	* @return true, wenn Ableitung benutzt werden soll
	*/
	virtual bool obj_diff(DiffStructure &s);

	virtual void integral(double *f, double t, const double *x, const double *u, const double *p);
	virtual bool integral_structure(DiffStructure &s);
	virtual bool integral_diff(DiffStructure &s, double t, const double *x, const double *u, const double *p);

	/** ODE-System / Dynamik */
	virtual void ode(double *dx, double t, const double *x, const double *u, const double *p)=0;
	/** Struktur des ODE-System */
	virtual bool ode_structure(DiffStructure &s);
	/** Ableitung des ODE-System */
	virtual bool ode_diff(DiffStructure &s, double t, const double *x, const double *u, const double *p);
	/** Ableitung des ODE-System nach freien Parameter */
	virtual bool ode_diff_p(DiffStructure &s, double t, const double *x, const double *u, const double *p, int index);

	/** Box-Beschraenkungen fuer Zustand */
	virtual void x_boundary(double *x_low, double *x_upp);
	/** Box-Beschraenkungen fuer Steuerung */
	virtual void u_boundary(double *u_low, double *u_upp);
	/** Box-Beschraenkungen fuer freie Parameter */
	virtual void p_boundary(double *p_low, double *p_upp);
	/** allgemeine Box-Beschraenkungen */
	virtual void var_boundary(double *x_low, double *x_upp);

	/** Randbedingungen */
	virtual void rand(double *r);
	/** Schranken fuer Randbedingungen */
	virtual void rand_boundary(double *r_low, double *r_upp);
	/** Struktur der Randbedingungen */
	virtual bool rand_structure(DiffStructure &s);
	/** Ableitung der Randbedingungen */
	virtual bool rand_diff(DiffStructure &s);

	/** Nebenbedingungen */
	virtual void neben(double *c, double t, const double *x, const double *u, const double *p);
	/** Schranken fuer Nebenbedingungen */
	virtual void neben_boundary(double *c_low, double *c_upp);
	/** Struktur der Nebenbedingungen */
	virtual bool neben_structure(DiffStructure &s);
	/** Ableitung der Nebenbedingungen */
	virtual bool neben_diff(DiffStructure &s, double t, const double *x, const double *u, const double *p);
	/** Ableitung der Nebenbedingungen nach freien Parametern */
	virtual bool neben_diff_p(DiffStructure &s, double t, const double *x, const double *u, const double *p, int index);

	/* --------------- Startschaetzung --------------- */
	/** wird vor der Optimierung einmalig aufgerufen */
	virtual void init();

	/** Initialisierung der Zustaende
	* @param x Zustand
	* @param i dis Punkt
	* @param dis Anz aller Gitterpunkte
	*/
	virtual void x_init(double *x, int i, int dis);
	/** Initialisierung der Steuerungen
	* @param u Steuerung
	* @param i dis Punkt
	* @param dis Anz aller Gitterpunkte
	*/
	virtual void u_init(double *u, int i, int dis);
	/** Initialisierung der freien Parameter
	* @param p freie Parameter
	*/
	virtual void p_init(double *p);

	virtual void zen_init(double *zen);

	/* --------------- Text-Ausgabe --------------- */
	/** wird nach Optimierung aufgerufen */
	virtual void terminate();
	/** wird nach jeden Optimierungsschritt aufgerufen */
	virtual bool step();

	/** Grafische Ausgabe */
#ifdef TRANSWORHP_GRAPHICS
	virtual void OpenWindows(Viewer *gr);
	virtual std::string GetXTitle(int d);
	virtual std::string GetUTitle(int d);

	/** Methode zur Auswahl welche Fenster dargestellt werden sollen */
	virtual void selectWindows(Viewer *viewer);
#endif

	virtual void localinit();


public:

	const std::string id; /**< Problemname */
	const int unique_id;

	std::vector<int> multinodes;
	const std::vector<int> box_neben;

	int n_dis; /**< Anzahl Gitterpunkte */
	const int n_ode; /**< Anzahl Zustaende */
	const int n_ctrl; /**< Anzahl Steuerungen */
	const int n_param; /**< Anzahl freier Parameter */
	const int n_rand; /**< Anzahl Randbedingungen */
	const int n_neben; /**< Anzahl Nebenbedingungen */
	const int n_integral; /**< Anzahl Integralterm */
	int n_zen;

	bool freetime; /**< freie Endzeit */

	double t0; /**< Startzeitpunkt */
	double tf; /**< Endzeitpunkt */


	double x(int dis, int ode) const;
	double u(int dis, int ctrl) const;
	double p(int param) const;
	double x__(int dis, int ode) const;
	double u__(int dis, int ctrl) const;
	int x_index(int dis, int ode) const;
	int x_index__(int dis, int ode) const;
	int u_index(int dis, int ctrl) const;
	int u_index__(int dis, int ctrl) const;
	int p_index(int param) const;
	int x_indexode(int ode) const;
	int u_indexode(int ctrl) const;
	int p_indexode(int param) const;
};




class DllExport ScaleTransWorhpProblem : public TransWorhpProblem {

public:

	ScaleTransWorhpProblem(const TWdimension &TWdata);

	void ode(double *dx, double t, const double *x, const double *u, const double *p) override final;
	bool ode_diff(DiffStructure &s, double t, const double *x, const double *u, const double *p) override final;
	bool ode_diff_p(DiffStructure &s, double t, const double *x, const double *u, const double *p, int index) override final;

	void x_init(double *x, int i, int dis) override final;
	void u_init(double *u, int i, int dis) override final;

	void x_boundary(double *x_low, double *x_upp) override final;
	void u_boundary(double *u_low, double *u_upp) override final;

	virtual void SC_ode(double *dx, double t, const double *x, const double *u, const double *p)=0;
	virtual bool SC_ode_diff(DiffStructure &s, double t, const double *x, const double *u, const double *p);
	virtual bool SC_ode_diff_p(DiffStructure &s, double t, const double *x, const double *u, const double *p, int index);

	virtual void SC_x_boundary(double *x_low, double *x_upp);
	virtual void SC_u_boundary(double *u_low, double *u_upp);

	virtual void SC_x_init(double *x, int i, int dis);
	virtual void SC_u_init(double *u, int i, int dis);

	virtual void u_scale(double *u);
	virtual void u_unscale(double *u);
	virtual void x_scale(double *x);
	virtual void x_unscale(double *x);
};

class DllExport EmptyTransWorhpProblem : public TransWorhpProblem {
public:
	EmptyTransWorhpProblem(const TWdimension &TWdata);

	double obj() override final;
	void ode(double *dx, double t, const double *x, const double *u, const double *p) override final;
};

}
