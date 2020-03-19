#pragma once

namespace tw {

/** Abstrakte Basisklasse fuer Splines */
class TWbaseSpline {
public:
	/** Destruktor */
	virtual ~TWbaseSpline(){}
  
	virtual double eval(double x)=0;
	virtual double deriv(double x)=0;
};

/** Klasse zum erstellen von linearen und kubischen "Splines"
* @author Matthias Rick
*/
class TWspline : public TWbaseSpline {
public:
	/** erstellt kubischen Spline in Intervall [x[k], x[k+1]] der Form: 
	* s(x) = a(x-xk)^3 + b(x-xk)^2 + c(x-xk) + d
	* @param xk Intervall Anfang
	* @param xk1 Intervall Ende
	* @param fxk Funktionswert an Intervall Anfang
	* @param fxk1 Funktionswert an Intervall Ende
	* @param dfxk Wert der Ableitung an Intervall Anfang
	* @param dfxk1 Wert der Ableitung an Intervall Ende
	*/
	TWspline(double xk, double xk1, double fxk, double fxk1, double dfxk, double dfxk1);
	/** erstellt linearen Spline in Intervall [x[k], x[k+1]] der Form:
	* s(x) = c(x-xk) + d
	* @param xk Intervall Anfang
	* @param xk1 Intervall Ende
	* @param fxk Funktionswert an Intervall Anfang
	* @param fxk1 Funktionswert an Intervall Ende
	*/
	TWspline(double xk, double xk1, double fxk, double fxk1);
	/** Destruktor */
	~TWspline();
	
	/** bestimmt erneut die Koeffizienten */
	void setCoef(double xk, double xk1, double fxk, double fxk1, double dfxk, double dfxk1);
	void setCoef(double xk, double xk1, double fxk, double fxk1);

	/** wertet den Spline an der Stelle x aus
	* @param x Stelle an der der Spline ausgewertet wird
	*/
	double eval(double x);
	/** wertet die Ableitung des Spline an der Stelle x aus
	* @param x Stelle an der die Ableitung ausgewertet wird
	*/
	double deriv(double x);

private:
	/** Koeffizienten */
	double a,b,c,d;

	/** Entwicklungspunkt */
	double xk;
};

/** Klasse zum erstellen von "Splines" 5. und 2. Grades fuer Hermite-Simpson-Verfahren
* @author Matthias Rick
*/
class TWsplineHS : public TWbaseSpline {
public:
	/** erstellt Spline in Intervall [x[k], x[k+1]] der Form: 
	* s(x) = a(x-xk)^5 + b(x-xk)^4 + c(x-xk)^3 + d(x-xk)^2 + e(x-xk) + f
	*
	* @param xk Intervall Anfang
	* @param xk05 Intervall Mitte
	* @param xk1 Intervall Ende
	*
	* @param fxk Funktionswert an Intervall Anfang
	* @param fxk05 Funktionswert an Intervall Mitte
	* @param fxk1 Funktionswert an Intervall Ende
	*
	* @param dfxk Wert der Ableitung an Intervall Anfang
	* @param dfxk05 Wert der Ableitung an Intervall Mitte
	* @param dfxk1 Wert der Ableitung an Intervall Ende
	*/
	TWsplineHS(double xk, double xk05, double xk1, double fxk, double fxk05, double fxk1, double dfxk, double dfxk05, double dfxk1);
	
	/** erstellt Spline in Intervall [x[k], x[k+1]] der Form:
	* s(x) = d(x-xk)^2 + e(x-xk) + f
	*
	* @param xk Intervall Anfang
	* @param xk05 Intervall Mitte
	* @param xk1 Intervall Ende
	*
	* @param fxk Funktionswert an Intervall Anfang
	* @param fxk05 Funktionswert an Intervall Mitte
	* @param fxk1 Funktionswert an Intervall Ende
	*/
	TWsplineHS(double xk, double xk05, double xk1, double fxk, double fxk05, double fxk1);
	
	/** Destruktor */
	~TWsplineHS();

	/** wertet den Spline an der Stelle x aus
	* @param x Stelle an der der Spline ausgewertet wird
	*/
	double eval(double x);
	/** wertet die Ableitung des Spline an der Stelle x aus
	* @param x Stelle an der die Ableitung ausgewertet wird
	*/
	double deriv(double x);

private:
	/** Koeffizienten */
	double a,b,c,d,e,f;

	/** Entwicklungspunkt */
	double xk;
};

}
