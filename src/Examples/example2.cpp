#include "TransWORHP.h"
#include "../Models/Model1v1.h"
#include "../Models/Model2v1.h"
#include <vector>
using namespace std;

int main() {

	double xo  = 100.;
	double yo  = 200.;
	double psi = 0;
	double u   = 7.75;
	double v   = 0;
	double r   = 0;

	double delta = 0.2;
	double n     = 1.5;

	double Uc    = 7.75;
	double psic  = M_PI;
	double Uw    = 7.75;
	double psiw  = -M_PI/2;
	double h     = 100.;

	double dx[6];

	double t = 0;

	double x[6] = {xo, yo, psi, u, v, r};

	double ctrl[2] = {delta, n};

	double q[5] = {Uc, psic, Uw, psiw, h};

	vector <BaseModel*> ship ( 2 );
	ship[0] = new Model1v1 ( "ParamXML/s1v1_m1v1_p1.xml" );
	ship[1] = new Model1v1 ( "ParamXML/s1v1_m1v1_p1.xml" );

	ship[0]->ode ( dx, x, ctrl, NULL, q );

	cout << endl << endl;
	cout << "Ship 1.1" << endl;
	cout << "Model 1.1" << endl;
	cout << "Parameter set 1" << endl;
	cout << endl;
	cout << setprecision ( 6 ) << fixed;
	cout << "Time derivative of xo : " << setw ( 10 ) << right << dx[0] << " m/s"     << endl;
	cout << "Time derivative of yo : " << setw ( 10 ) << right << dx[1] << " m/s"     << endl;
	cout << "Time derivative of psi: " << setw ( 10 ) << right << dx[2] << " rad/s"   << endl;
	cout << "Time derivative of u  : " << setw ( 10 ) << right << dx[3] << " m/s^2"   << endl;
	cout << "Time derivative of v  : " << setw ( 10 ) << right << dx[4] << " m/s^2"   << endl;
	cout << "Time derivative of r  : " << setw ( 10 ) << right << dx[5] << " rad/s^2" << endl;
	cout << endl << endl;


	ship[1]->ode ( dx, x, ctrl, NULL, q );

	cout << endl << endl;
	cout << "Ship 1.1" << endl;
	cout << "Model 2.1" << endl;
	cout << "Parameter set 1" << endl;
	cout << endl;
	cout << setprecision ( 6 ) << fixed;
	cout << "Time derivative of xo : " << setw ( 10 ) << right << dx[0] << " m/s"     << endl;
	cout << "Time derivative of yo : " << setw ( 10 ) << right << dx[1] << " m/s"     << endl;
	cout << "Time derivative of psi: " << setw ( 10 ) << right << dx[2] << " rad/s"   << endl;
	cout << "Time derivative of u  : " << setw ( 10 ) << right << dx[3] << " m/s^2"   << endl;
	cout << "Time derivative of v  : " << setw ( 10 ) << right << dx[4] << " m/s^2"   << endl;
	cout << "Time derivative of r  : " << setw ( 10 ) << right << dx[5] << " rad/s^2" << endl;
	cout << endl << endl;


	return 0;
}





















