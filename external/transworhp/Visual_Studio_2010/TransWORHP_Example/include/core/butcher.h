#ifndef butcher_h
#define butcher_h

#include "../base/defines.h"
#include <iostream>

#ifdef WIN32
#ifndef MINGW
#define DllExport __declspec( dllexport )
#pragma warning (disable : 4251)
#else
#define DllExport
#endif
#else 
#define DllExport
#endif


class TransWorhp;
class ExplTransWorhp;

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
	Butcher(int n_ode, int n_ctrl);
	~Butcher();
	
	/** Einlesen eines vordefinierten Tableaus. */
	//void Init(int index, double stepsize=1e-3, double abserr=1e-6, double relerr=1e-6);
	void Init(int index, double stepsize=1e-3);

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
	int RungeKutta(TransWorhp *tw, double &t, const double t0, const double tend, double *x, const double *u1, const double *u2, const double *param, double &h);
	/** externes RK-Verfahren */
	int RungeKuttaF(TransWorhp *tw, double &t, double t0, double tend, double *x, double *u1, double *u2, double *param, int startflag);

	bool hasStepSize() const {return p!=0;}
	int getStepSize() const {return p;}
	void setStepSize(int i) {p = i;}
	int stufen() const {return s;}

	/** Bezeichnungen der Butchertableaus */
	static const std::string names[20];
	/** Butchertableaus */
	static const std::string tableaus[20];
	
	double h0; /**< Schrittweite für Verfahren ohne SW-Steuerung */

	double abserr;
	double relerr;

private:
	// speichern des Tableaus
	double a[20][20];
	double c[20];
	double b[20];
	double b2[20];
	
	double *k[20];
	
	const double amin;
	const double amax;
	const double aa;
	
	double maxEpsEst;
	double minEpsEst;
	
	const double eps;
	const double hmin;
	
	
	double *tmp_x, *tmp_u;
	double *tmp_ret1, *tmp_ret2;

	/** Stufen */
	int s;
	/** Ordnung des Verfahrens. wenn 0: Verfahren ohne SW-Steuerung */
	int p;
	
	/** Interpretatation von rationalen Zahlen der Form "2/3" oder "2". */
	double getratio(std::istream &os);
};


#endif