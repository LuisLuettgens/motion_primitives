#include "TWspline.h"

#include <iostream>
#include <cmath>

namespace tw {

TWspline::TWspline(double xk_, double xk1, double fxk, double fxk1, double dfxk, double dfxk1) {
	
	setCoef(xk_, xk1, fxk, fxk1,dfxk, dfxk1);
}

TWspline::TWspline(double xk_, double xk1, double fxk, double fxk1) {
	
	setCoef(xk_, xk1, fxk, fxk1);
}

TWspline::~TWspline() {}

void TWspline::setCoef (double xk_, double xk1, double fxk, double fxk1, double dfxk, double dfxk1 ) {
	
	xk = xk_;
	
	const double h = xk1 - xk;
	
	a = ( dfxk1*h + dfxk*h + 2*fxk - 2*fxk1 ) / pow(h, 3);
	//a = (dfxk1*h + dfxk*h + 2 * fxk - 2 * fxk1) / (h*h*h);
	b = ( dfxk1 - 3 * a*pow(h, 2) - dfxk ) / (2 * h);
	//b = -( h*(2*dfxk+dfxk1) + 3*(fxk-fxk1) ) / (pow(h,2));
	//b = -(h*(2 * dfxk + dfxk1) + 3 * (fxk - fxk1)) / (h*h);
	c = dfxk;
	d = fxk;
}

void TWspline::setCoef (double xk_, double xk1, double fxk, double fxk1 ) {
	
	xk = xk_;
	
	a = 0;
	b = 0;
	c = (fxk1 - fxk) / (xk1 - xk);
	d = fxk;
}

double TWspline::eval(double x) {
	const double h = x - xk;
	return a*h*h*h + b*h*h + c*h + d;
}

double TWspline::deriv(double x) {
	const double h = x - xk;
	return 3*a*h*h + 2*b*h + c;
}


////////////////////////////////////////////////////////////////////

/* xk05 wird nicht benoetigt, da xk05=(xk+xk1)/2 */
TWsplineHS::TWsplineHS(double xk_, double /*xk05*/, double xk1, double fxk, double fxk05, double fxk1, double dfxk, double dfxk05, double dfxk1) : xk(xk_) {
	const double h = xk1 - xk;
	
	// mit Maple
	a = (4*(4*h*dfxk05-6*fxk1+h*dfxk1+h*dfxk+6*fxk))/(h*h*h*h*h);
	b = -(4*(-13*fxk1+2*h*dfxk1+3*h*dfxk+17*fxk-4*fxk05+10*h*dfxk05))/(h*h*h*h);
	c = (5*h*dfxk1+13*h*dfxk+66*fxk-32*fxk05+32*h*dfxk05-34*fxk1)/(h*h*h);
	d = -(8*h*dfxk05-7*fxk1+h*dfxk1+6*h*dfxk+23*fxk-16*fxk05)/(h*h);
	e = dfxk;
	f = fxk;
	//std::cout << "spline5HS: " << a << ", " << b << ", " << c << ", " << d << ", " << e << ", " << f << std::endl;
}

TWsplineHS::TWsplineHS(double xk_, double /*xk05*/, double xk1, double fxk, double fxk05, double fxk1) : xk(xk_) {
	const double h = xk1 - xk;
  
	a = 0;
	b = 0;
	c = 0;
	d = (2*(fxk-2*fxk05+fxk1))/(h*h);
	e = -(3*fxk-4*fxk05+fxk1)/h;
	f = fxk;
	//std::cout << "spline2HS: " << a << ", " << b << ", " << c << ", " << d << ", " << e << ", " << f << std::endl;
}

TWsplineHS::~TWsplineHS() {}


double TWsplineHS::eval(double x) {
	const double h = x - xk;
	return a*h*h*h*h*h + b*h*h*h*h + c*h*h*h + d*h*h + e*h + f;
}

double TWsplineHS::deriv(double x) {
	const double h = x - xk;
	return 5*a*h*h*h*h + 4*b*h*h*h + 3*c*h*h + 2*d*h + e;
}

}
