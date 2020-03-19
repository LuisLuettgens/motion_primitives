/*----------------------------------------------------------------
 *
 *  Example: Emergency landing of a Hypersonic Flight System
 *
 *----------------------------------------------------------------*/

#ifdef WIN32
#include "windows.h"
#endif

#include "../src/core/TransWORHP.h"

#include <fstream>
#include <sstream>

using namespace std;

void ParseString(std::string& doubleList, std::vector<double>& vDoubleOut)
{
  std::stringstream strm;
  strm.str(doubleList);
  double curVal = 0;
  while(strm >> curVal)
    {
      vDoubleOut.push_back(curVal);
    }
}

class EmergencyLanding : public TransWorhp {
public:

	EmergencyLanding(int dis) : TransWorhp("Emergency landing", dis, 6, 2, 1, 0, 0) {
		freetime = 1;

		F = 305.0;
		r0 = 6.371e6;
		CD0 = 0.017;
		k = 2.0;
		rho0 = 1.249512; // Unperturbed case
		beta = 1.0 / 6900.0;
		g0 = 9.80665;
		omega = 7.27e-5;
		m = 115000.0;
		scaleA = 10000.; // Scaling for x[3] - altitude
		scaleV = 1000.; // scaling for x[0] - velocity
	}

	string GetXTitle(int d) override {
		if (d==0) return "v - velocity";
		if (d==1) return "gamma - inclination";
		if (d==2) return "chi - azimut";
		if (d==3) return "h - height";
		if (d==4) return "lambda - latitude";
		if (d==5) return "theta - longitude";
	}

	string GetUTitle(int d) override {
		if (d==0) return "C - drag";
		if (d==1) return "mu - roll angle";
	}

        void x_init(double *x, int i, int dis) override {
	  x[0] = 2150.54529 / scaleV;
	  x[1] = 0.152018177;
	  x[2] = 2.2689279889;//95.0*M_PI/180.0;
	  x[3] = 33900.0 / scaleA;
	  x[4] = 0.8651597102;
	  x[5] = 0.1980948701;
	}

	void p_init(double *p) override {
	  p[0] = 728.85;
	}

	void u_init(double *u, int i, int dis) override {
		double myAux = i/(dis-1.);
		u[1]	= (myAux - 1)*(myAux - 1) * 0.8;
		u[0] = 0.1;
	// u[1]  = 0.8 - 0.8*i/(dis-1);
	}

	double obj() override {
		// Objective function
		return -( (x(n_dis-1, 4) - x(0, 4) ))* ( (x(n_dis-1,4) - x(0, 4) ))
 - ( (x(n_dis-1, 5) - x(0, 5) )) * ( (x(n_dis-1,5) - x(0, 5) ));
	}

	bool obj_structure(DiffStructure &s) override {
	  s(0, x_index(n_dis-1, 4) );
	  s(0, x_index(n_dis-1, 5) );
	  return true;
	}

	bool obj_diff(DiffStructure &s) override {
	  s(0, x_index(n_dis-1, 4) ) = - 2.0 * ( (x(n_dis-1, 4) - x(0, 4) ));
	  s(0, x_index(n_dis-1, 5) ) = - 2.0 * ( (x(n_dis-1, 5) - x(0, 5) ));
	  return true;
	}

	void ode(double *dx, double t, const double *x, const double *u, const double *p) override {
		// Dependencies in comments
		double rho; // x[3]
		double q; // x[0] ,x[3]
		double L; // x[0], x[3], u[0]
		double CD; // u[0]
		double D; // x[0], x[3], u[0]
		double R; // x[3]
		double g; // x[3]
		CalculateAuxiliaries(x, u, rho, q, L, CD, D, R, g);

		// Some optimisation
		double sinx1 = sin(x[1]);
		double cosx1 = cos(x[1]);
		double tanx1 = sinx1/cosx1;
		double sinx2 = sin(x[2]);
		double cosx2 = cos(x[2]);
		double sinx4 = sin(x[4]);
		double cosx4 = cos(x[4]);
		double tanx4 = sinx4/cosx4;
		double sinu1 = sin(u[1]);
		double cosu1 = cos(u[1]);
		double x0 = x[0] * scaleV;
		double x3 = x[3] * scaleA;

		// Model
		dx[0] = -D/m - g*sinx1 + omega*omega * cosx4*(sinx1*cosx4 - cosx1*sinx2*sinx4)*R;
		dx[1] = L*cosu1/(m*(x0)) - (g/(x0) - (x0)/R)*cosx1 + 2.0*omega*cosx2*cosx4 + omega*omega * cosx4*(sinx1*sinx2*sinx4 + cosx1*cosx4)*R/(x0);
		dx[2] = L*sinu1/(m*(x0)*cosx1) - cosx1*cosx2*tanx4*(x0)/R + 2.0*omega*(sinx2*cosx4*tanx1-sinx4) - omega*omega*cosx4*sinx4*cosx2*R/((x0)*cosx1);
		dx[3] = (x0)*sinx1;
		dx[4] = cosx1*sinx2*(x0)/R;
		dx[5] = cosx1*cosx2*(x0)/(R*cosx4);

		// scaling
		dx[0] /= scaleV;
		dx[3] /= scaleA;

		// free final time
		dx[0] *= p[0];
		dx[1] *= p[0];
		dx[2] *= p[0];
		dx[3] *= p[0];
		dx[4] *= p[0];
		dx[5] *= p[0];
	}

	bool ode_structure(DiffStructure &s) override {
		s(0, x_indexode(0));
		s(0, x_indexode(1));
		s(0, x_indexode(2));
		s(0, x_indexode(3));
		s(0, x_indexode(4));
		s(0, u_indexode(0));

		s(1, x_indexode(0));
		s(1, x_indexode(1));
		s(1, x_indexode(2));
		s(1, x_indexode(3));
		s(1, x_indexode(4));
		s(1, u_indexode(0));
		s(1, u_indexode(1));

		s(2, x_indexode(0));
		s(2, x_indexode(1));
		s(2, x_indexode(2));
		s(2, x_indexode(3));
		s(2, x_indexode(4));
		s(2, u_indexode(0));
		s(2, u_indexode(1));

		s(3, x_indexode(0));
		s(3, x_indexode(1));

		s(4, x_indexode(0));
		s(4, x_indexode(1));
		s(4, x_indexode(2));
		s(4, x_indexode(3));

		s(5, x_indexode(0));
		s(5, x_indexode(1));
		s(5, x_indexode(2));
		s(5, x_indexode(3));
		s(5, x_indexode(4));

		// Derivative w.r.t. p[0]
		s(0, p_indexode(0) );
		s(1, p_indexode(0) );
		s(2, p_indexode(0) );
		s(3, p_indexode(0) );
		s(4, p_indexode(0) );
		s(5, p_indexode(0) );

		return true;
	}

	bool ode_diff(DiffStructure& s, double t, const double* x, const double* u, const double* p) override {
		// Double checked those, but TransWORHP prefers Finite Differences.
		// Something still buggy down here. Working fine with finite difference though.
		return false;
/*		double rho; // x[3]
		double q; // x[0] ,x[3]
		double L; // x[0], x[3], u[0]
		double CD; // u[0]
		double D; // x[0], x[3], u[0]
		double R; // x[3]
		double g; // x[3]
		CalculateAuxiliaries(x, u, rho, q, L, CD, D, R, g);

		// Some optimisation
		double sinx1 = sin(x[1]);
		double cosx1 = cos(x[1]);
		double tanx1 = sinx1/cosx1;
		double sinx2 = sin(x[2]);
		double cosx2 = cos(x[2]);
		double sinx4 = sin(x[4]);
		double cosx4 = cos(x[4]);
		double tanx4 = sinx4/cosx4;
		double sinu1 = sin(u[1]);
		double cosu1 = cos(u[1]);
		double x0 = x[0] * scaleV;
		double x3 = x[3] * scaleA;

		double dRdX3 = 1.;
		double drhodX3 = -beta * rho0 * exp(-beta*x3);
		double dqdX3 = 0.5 * drhodX3 * (x0) * (x0);
		double dDdX3 = dqdX3 * F * CD;
		double dgdX3 = -2. * g0 * r0 * r0 * (1. / ((r0 + x3) * (r0 + x3) * (r0 + x3)) ) * dRdX3;

		double dDdU0 = 0.5 * rho * x0 * x0 * F * (2. * k * u[0]);
		
		double dLdX3 = dqdX3 * F * u[0];
		double dLdU0 = q * F;

		s(0, x_indexode(0)) = -1./m * rho * (x0) * ( F * CD );
		s(0, x_indexode(1)) = - g * cosx1 + omega*omega * cosx4*(cosx1*cosx4 + sinx1*sinx2*sinx4) * R;
		s(0, x_indexode(2)) = omega * omega * cosx4 * ( - cosx1*cosx2*sinx4) * R;
		s(0, x_indexode(3)) = -dDdX3/m - dgdX3 * sinx1 + omega*omega * cosx4 * ( sinx1 * cosx4 - cosx1 * sinx2 * sinx4 ) * dRdX3;
		s(0, x_indexode(4)) = (omega * omega * R) * (- 2. * cosx4 * sinx4 * sinx1 - cosx1 * sinx2 * ( - sinx4 * sinx4 + cosx4 * cosx4));
		s(0, u_indexode(0)) = -1./m * dDdU0;

		s(1, x_indexode(0)) = 0.5 * rho0 * exp(-beta*x3) * F * u[0] * cosu1 / m - cosx1 * ( ( - g / ((x0) * (x0)) ) - 1. / ( R )) - omega*omega * cosx4*(sinx1*sinx2*sinx4 + cosx1*cosx4)*R/((x0)*(x0));
		s(1, x_indexode(1)) = (g/(x0) - (x0)/R)*sinx1  + omega*omega * cosx4*(cosx1*sinx2*sinx4 - sinx1*cosx4)*R/(x0);
		s(1, x_indexode(2)) =  - 2.0*omega*sinx2*cosx4 + omega*omega * cosx4*(sinx1*cosx2*sinx4 )*R/(x0);
		s(1, x_indexode(3)) = dLdX3 * cosu1/(m*(x0)) - cosx1 * (dgdX3 / (x0) + ( (x0) / (R * R) ) * dRdX3 )  + omega*omega * cosx4*(sinx1*sinx2*sinx4 + cosx1*cosx4)/(x0);
		s(1, x_indexode(4)) = - 2.0 * omega * cosx2 * sinx4 + omega * omega * ((sinx1 * sinx2 * (-sinx4*sinx4 + cosx4 * cosx4)) + cosx1 * (-2*cosx4*sinx4)) * R/(x0);
		s(1, u_indexode(0)) = dLdU0 * cosu1/(m*(x0));
		s(1, u_indexode(1)) = - L * sinu1 / (m*(x0));

		s(2, x_indexode(0)) = 0.5 * rho0 * exp(-beta*x[3]) * F * u[0] * sinu1 / (m * cosx1) - cosx1*cosx2*tanx4/R + omega*omega*cosx4*sinx4*cosx2*R/((x0)*(x0)*cosx1);
		s(2, x_indexode(1)) = L*sinu1*sinx1/(m*(x0)*cosx1*cosx1) + sinx1*cosx2*tanx4*(x0)/R + 2.0*omega*(sinx2*cosx4/(cosx1*cosx1)) - omega*omega*cosx4*sinx4*cosx2*R*sinx1/((x0)*cosx1*cosx1);
		s(2, x_indexode(2)) = cosx1*sinx2*tanx4*(x0)/R + 2.0*omega*(cosx2*cosx4*tanx1) + omega*omega*cosx4*sinx4*sinx2*R/((x0)*cosx1);
		s(2, x_indexode(3)) = dLdX3*sinu1/(m*(x0)*cosx1) + cosx1*cosx2*tanx4*(x0)/(R*R)  - omega*omega*cosx4*sinx4*cosx2/((x0)*cosx1);
		s(2, x_indexode(4)) = - cosx1*cosx2*(1./(cosx4*cosx4))*(x0)/R + 2.0*omega*(sinx2*(-1.)*sinx4*tanx1-cosx4)  - omega*omega*(cosx4*cosx4*(-sinx4 * sinx4))*cosx2*R/((x0)*cosx1);
		s(2, u_indexode(0)) = dLdU0 * sinu1/(m*(x0)*cosx1);
		s(2, u_indexode(1)) = L * cosu1/(m*(x0)*cosx1);

		s(3, x_indexode(0)) = ( sinx1 );
		s(3, x_indexode(1)) = ( (x0)*cosx1 );


		s(4, x_indexode(0)) = cosx1*sinx2/R;
		s(4, x_indexode(1)) = -sinx1*sinx2*(x0)/R;
		s(4, x_indexode(2)) = cosx1*cosx2*(x0)/R;
		s(4, x_indexode(3)) = -cosx1*sinx2*(x0)/(R*R);

		s(5, x_indexode(0)) = cosx1*cosx2/(R*cosx4);
		s(5, x_indexode(1)) = -sinx1*cosx2*(x0)/(R*cosx4);
		s(5, x_indexode(2)) = -cosx1*sinx2*(x0)/(R*cosx4);
		s(5, x_indexode(3)) = -cosx1*cosx2*(x0)/(R*R*cosx4);
		s(5, x_indexode(4)) = cosx1*cosx2*(x0)*sinx4/(R*cosx4*cosx4);

		//scaling
		s(0, x_indexode(0)) /= scaleV;
		s(0, x_indexode(1)) /= scaleV;
		s(0, x_indexode(2)) /= scaleV;
		s(0, x_indexode(3)) /= scaleV;
		s(0, x_indexode(4)) /= scaleV;
		s(0, u_indexode(0)) /= scaleV;
		s(3, x_indexode(0)) /= scaleA;
		s(3, x_indexode(1)) /= scaleA;
		return true; */
	}

	bool ode_diff_p(DiffStructure& s, double t, const double *x, const double *u, const double *p, int index) override {
		// Dependencies in comments
		double rho; // x[3]
		double q; // x[0] ,x[3]gg107
		double L; // x[0], x[3], u[0]
		double CD; // u[0]
		double D; // x[0], x[3], u[0]
		double R; // x[3]
		double g; // x[3]
		CalculateAuxiliaries(x, u, rho, q, L, CD, D, R, g);

		// Some optimisation
		double sinx1 = sin(x[1]);
		double cosx1 = cos(x[1]);
		double tanx1 = sinx1/cosx1;
		double sinx2 = sin(x[2]);
		double cosx2 = cos(x[2]);
		double sinx4 = sin(x[4]);
		double cosx4 = cos(x[4]);
		double tanx4 = sinx4/cosx4;
		double sinu1 = sin(u[1]);
		double cosu1 = cos(u[1]);
		double x0 = x[0] * scaleV;
		double x3 = x[3] * scaleA;

		s(0, p_indexode(0) ) = -D/m - g*sinx1 + omega*omega * cosx4*(sinx1*cosx4 - cosx1*sinx2*sinx4)*R;
		s(1, p_indexode(0) ) = L*cosu1/(m*(x0)) - (g/(x0) - (x0)/R)*cosx1 + 2.0*omega*cosx2*cosx4 + omega*omega * cosx4*(sinx1*sinx2*sinx4 + cosx1*cosx4)*R/(x0);
		s(2, p_indexode(0) ) = L*sinu1/(m*(x0)*cosx1) - cosx1*cosx2*tanx4*(x0)/R + 2.0*omega*(sinx2*cosx4*tanx1-sinx4) - omega*omega*cosx4*sinx4*cosx2*R/((x0)*cosx1);
		s(3, p_indexode(0) ) = (x0)*sinx1;
		s(4, p_indexode(0) ) = cosx1*sinx2*(x0)/R;
		s(5, p_indexode(0) ) = cosx1*sinx2*(x0)/(R*cosx1);

		// scaling
		s(0, p_indexode(0) ) /= scaleV;
		s(3, p_indexode(0) ) /= scaleA;
		return true;
	}

	void x_boundary(double *x_low, double *x_upp) override {
		x_low[0] = 100 / scaleV; // minimum velocity
		x_upp[0] = 4000 / scaleV;
		x_low[1] = -M_PI; // angles
		x_upp[1] = M_PI;
		x_low[2] = -M_PI;
		x_upp[2] = M_PI;
		x_low[3] = 500 / scaleA; // altitude
		x_upp[3] = 100000 / scaleA;
		x_low[4] = -1.0; // latitude
		x_upp[4] = 1.0;
		x_low[5] = -1.0; // longitude
		x_upp[5] = 1.0;
	}

	void u_boundary(double *u_low, double *u_upp) override {
		u_low[0] = 0.01; // Bounds from PDF
		u_upp[0] = 0.18326;
		u_low[1] = -M_PI/2.0;
		u_upp[1] = M_PI/2.0;

	}

	void p_boundary(double *p_low, double *p_upp) override {
		p_low[0] = 700; // Useful bounds to time
		p_upp[0] = 800;
	}

	void var_boundary(double *x_low, double *x_upp) override {

		// Initial boundary conditions
		x_low[x_index(0,0)] = x_upp[x_index(0,0)] = 2150.54529 / scaleV;
		x_low[x_index(0,1)] = x_upp[x_index(0,1)] = 0.152018177;
		x_low[x_index(0,2)] = x_upp[x_index(0,2)] = 2.2689279889;//95.0*M_PI/180.0;
		x_low[x_index(0,3)] = x_upp[x_index(0,3)] = 33900.0 / scaleA;
		x_low[x_index(0,4)] = x_upp[x_index(0,4)] = 0.8651416299;
		x_low[x_index(0,5)] = x_upp[x_index(0,5)] = 0.1980907302;

		// Final boundary condition
		x_low[x_index(n_dis-1,3)] = x_upp[x_index(n_dis-1,3)] = 500.0 / scaleA;
	}
	
private:
	/* model constants */
	double F;
	double r0;
	double CD0;
	double k;
	double rho0; // Unperturbed case
	double beta;
	double g0;
	double omega;
	double m;
	
	double scaleA;
	double scaleV;

	void CalculateAuxiliaries(const double * x, const double * u, double & rho, double & q, double & L, double & CD, double & D, double & R, double & g) {
		double x0 = x[0] * scaleV;
		double x3 = x[3] * scaleA;

		rho = rho0 * exp(-beta*x3); // x[3]
		q = 0.5 * rho * (x0) * (x0); // x[0] ,x[3]
		L = q * F * u[0]; // x[0], x[3], u[0]
		CD = CD0 + k*u[0]*u[0]; // u[0]
		D = q * F * CD; // x[0], x[3], u[0]
		R = r0 + x3; // x[3]
		g = g0 * (r0 / R) * (r0 / R); // x[3]
	}
};

int main(int argv, char* argc[]) {

	TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv,argc);
	TWfolder folder(&twparameter,0);
	EmergencyLanding ph(twparameter.NDIS);
	folder.Add(&ph);

	Viewer *viewer = nullptr;
	if (twparameter.PLOT) viewer = new Viewer(&twparameter);

	folder.Init();
	folder.Init(viewer);


	// Open file containing good initial guess
/*	ifstream f;
	string s;
	vector<double> myData;
	f.open("emerLandInitGuess.dat", ios::in);
	getline(f, s); // Eat headline
	getline(f, s); // twice
	while (!f.eof()) {
	  getline(f, s);
	  ParseString(s, myData);
	}
	f.close();

	// First 81 entries correspond to old time parameter in [0, 1]
	if (ph.n_dis == 81) {
		for (size_t i = 0; i < ph.n_dis; ++i) {
			for (size_t j = 0; j < 6; ++j) {
				ph.X[ph.x_index(i, j)] = myData[ph.n_dis*(j+1) + i];
			}
		}
	}
	ph.ToMATLAB("test.m"); */
	//ph.FromMATLAB("test.m");
	ph.Integrate(2);
	folder.meshRef();

	delete viewer;

	return 0;
}
