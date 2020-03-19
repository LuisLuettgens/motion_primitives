#pragma once

#include "../base/defines.h"
#include "diffstructure.h"

#include <vector>
#include <iostream>

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

#if defined(TW_WITH_SUPERLU) || defined(WITH_LAPACK)
class ExplTransWorhp;
#endif


/** Butchertableau zur ODE-Integration.
 *
 *  Hier: Fest in names[] und tableaus[] codierte Verfahren
 *
 * Verfahren        |  s  |  p  |  q  |
 * HeunEuler        |  2  |  2  |  1  |
 * BogackiShampine  |  4  |  3  |  2  |
 * Uebung           |  4  |  2  |  3  |
 * Fehlberg         |  6  |  4  |  5  |
 * CashKarp         |  6  |  5  |  4  |
 * DormandPrince    |  7  |  4  |  5  |
 */
class DllExport Butcher {
public:

	/** Konstruktor. */
	Butcher(int n_ode, int n_ctrl, double abserr, double relerr);

	/** Einlesen eines vordefinierten Tableaus.
	@param index Butcher-Tableaus
	@param stepsize Schrittweite
	@param verbose 1=mit Ausgabe, 0=ohne Ausgabe
	*/
	void Init(int index, double stepsize=1e-3, bool verbose=true);
	//void Init(int index, double stepsize=1e-3, double abserr=1e-6, double relerr=1e-6);

	/** Einlesen eines Tableaus. */
	friend std::istream &operator>>(std::istream &os, Butcher &b);

	/** Ausgabe eines Tableaus. */
	friend std::ostream &operator<<(std::ostream &os, const Butcher &b);

	/** Ein Schritt der expliziten Integration.
	*
	* @param tw Pointer auf TW
	* @param t aktueller Zeitpunkt
	* @param t0 Intervallanfang
	* @param tend Intervallende
	* @param x Zustaende
	* @param u1 Steuerung am linken Rand
	* @param u2 Steuerung am rechten Rand
	* @param param Parametervektor
	* @param h aktuelle Schrittweite
	* @return -1: Schrittweite zu klein, 0: erfolgreicher Schritt, 1: Fehler fuer Schrittweite zu gross (erneuter Schritt noetig)
	*/
	int RungeKutta(TransWorhpProblem *tw, double &t, double t0, double tend, double *x, const double *u1, const double *u2, const double *param, double &h);
	/** externes RK-Verfahren */
	int RungeKuttaF(TransWorhpProblem *tw, double &t, double t0, double tend, double *x, double *u1, double *u2, double *param, int startflag);
	
#ifdef TW_WITH_BOOST
	int integrate_with_boost(TransWorhpProblem* tw, double t0, double tend, std::vector<double> &x, const double* param, int butchertableau) const;
#endif

#if defined(TW_WITH_SUPERLU) || defined(WITH_LAPACK)
	int ROW(ExplTransWorhp *tw, const double t0, const double tend, double *x, const double *u1, const double *u2, const double *param);
#endif

#ifdef TW_WITH_SUPERLU
	/** linear-implizites Rosenbrock-Wanner-Verfahren */
	int ROW_SUPERLU(ExplTransWorhp *tw, const double t0, const double tend, double *x, const double *u1, const double *u2, const double *param);
#endif

#ifdef WITH_LAPACK
	/** linear-implizites Rosenbrock-Wanner-Verfahren */
	int ROW_LAPACK(ExplTransWorhp *tw, const double t0, const double tend, double *x, const double *u1, const double *u2, const double *param);
#endif

	bool hasStepSize() const {return p!=0;}
	int getStepSize() const {return p;}
	void setStepSize(int i) {p = i;}
	int stufen() const {return s;}

	/** Bezeichnungen der Butchertableaus */
	static const std::vector<std::string> names;
	/** Butchertableaus */
	static const std::vector<std::string> tableaus;

	double h0; /**< Schrittweite für Verfahren ohne SW-Steuerung */

private:

	bool isInit;

	// speichern des Tableaus
	double a[20][20];
	double c[20];
	double b[20];
	double b2[20];

#if defined(TW_WITH_SUPERLU) || defined(WITH_LAPACK)
	//fuer ROW-Verfahren
	double gamma[2][2];
	double d[2];

	int nnzOfA;
#endif

	std::vector<std::vector<double>> k;

	const double amin;
	const double amax;
	const double aa;

	double maxEpsEst;
	double minEpsEst;

	const double eps;
	const double hmin;

	double abserr;
	double relerr;

#if defined(TW_WITH_SUPERLU) || defined(WITH_LAPACK)
	DiffStructure ROW_DS;
	bool ROW_DS_check;
#endif

	std::vector<double> tmp_x, tmp_u;
	std::vector<double> tmp_ret1, tmp_ret2;

	/** Stufen */
	int s;
	/** Ordnung des Verfahrens. wenn 0: Verfahren ohne SW-Steuerung */
	int p;

	const int n_ode;
	const int n_ctrl;

#ifdef TW_WITH_BOOST
	template<class stepper_type>
	int integrate_with_boost(TransWorhpProblem* tw, double t0, double tend, std::vector<double> &x, const double* param) const;

	template<class stepper_type>
	int integrate_with_boost_controlled(TransWorhpProblem* tw, double t0, double tend, std::vector<double> &x, const double* param) const;
#endif

	/** Interpretatation von rationalen Zahlen der Form "2/3" oder "2". */
	double getratio(std::istream &os);
};

}
