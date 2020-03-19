/*----------------------------------------------------------------
 *
 *  Example: Notlande
 *
 *----------------------------------------------------------------*/

#ifdef WIN32
#include "windows.h"
#endif

#include "TransWORHP.h"

using namespace std;


class NotlandePhase : public TransWorhp {
public:
	

	NotlandePhase(int dis) : TransWorhp("Notlande", dis,6,2,1,0,0) {
		freetime = 1;

	}


	string GetXTitle(int d) {
		if (d==0) return "v - velocity";
		if (d==1) return "gamma - inclination";
		if (d==2) return "chi - azimut";
		if (d==3) return "h - height";
		if (d==4) return "lambda - latitude";
		if (d==5) return "theta - longitude";
	}

	string GetUTitle(int d) {
		if (d==0) return "C - drag";
		if (d==1) return "mu - roll angle";
	}
	
	
	void p_init(double *p) {
		//p[0] = 700;
	}

	void u_init(double *u) {
		//u[0] = 0.1;
		//u[1] = 0.8;
	}
	
	double obj() {

		// Zielfunktion
		return -((x(n_dis-1,4) - x(0, 4))/x(0, 4))*((x(n_dis-1,4) - x(0, 4))/x(0, 4)) - ((x(n_dis-1,5) - x(0, 5))/x(0, 5))*((x(n_dis-1,5) - x(0, 5))/x(0, 5));
	}

	bool obj_structure(DiffStructure &s) {
		//s(0,x_index(n_dis-1,5));
		//s(0,p_index(0));
		return false;
	}
	
	bool obj_diff(DiffStructure &s) {
		//s(0,x_index(n_dis-1,5)) = 0.0;
		//s(0,p_index(0)) = 1;
		return false;
	}
	

	void ode(double *dx, double t, const double *x, const double *u, const double *p) {
		

		// PHYSIKALISCHE GROESSEN 
/*
		double R     = R0+x[3];
		double GH    = G0*R0*R0/(R*R);
		double RHOH  = RHO0*exp(-BETA*x[3]);
		double QH    = 0.5*RHOH*x[0]*x[0];
		double CD    = CD0+K*u[0]*u[0];

		double HILF1 = QH*FF;
		double D     = HILF1*CD;
		double L     = HILF1*u[0];
*/

		

		double D = 0.5*(1.249512 + 1.249512*1/100.0)*exp(-1.0/6900.0*x[3])*x[0]*x[0]*305.0*0.017+2*u[0]*u[0];
		double L = 0.5*(1.249512 + 1.249512*1/100.0)*exp(-1.0/6900.0*x[3])*x[0]*x[0]*305.0*u[0];
		double R = 6.371e+6+x[3];
		double g = 9.81*6.371e+6*6.371e+6/(R*R);
		double omega = 7.27e-5;
		double m = 115000;
		
		// Modell 
		dx[0] = -D/m - g*sin(x[1]) + omega*omega * cos(x[4])*(sin(x[1])*cos(x[4]) - cos(x[1])*sin(x[2])*sin(x[4]))*R;
		dx[1] = L*cos(u[1])/(m*x[0]) - (g/x[0] - x[0]/R)*cos(x[1]) + 2.0*omega*cos(x[2])*cos(x[4]) + omega*omega * cos(x[4])*(sin(x[1])*sin(x[2])*sin(x[4]) + cos(x[1])*cos(x[4]))*R/x[0];
		dx[2] = L*sin(u[1])/(m*x[0]*cos(x[1])) - cos(x[1])*cos(x[2])*tan(x[4])*x[0]/R + 2.0*omega*(sin(x[2])*cos(x[4])*tan(x[1])-sin(x[4])) - omega*omega*cos(x[4])*sin(x[4])*cos(x[2])*R/(x[0]*cos(x[1]));
		dx[3] = x[0]*sin(x[1]);
		dx[4] = cos(x[1])*sin(x[2])*x[0]/R;
		dx[5] = cos(x[1])*sin(x[2])*x[0]/(R*cos(x[1]));

		// freie Endzeit
		dx[0] = dx[0]*p[0];
		dx[1] = dx[1]*p[0];
		dx[2] = dx[2]*p[0];
		dx[3] = dx[3]*p[0];
		dx[4] = dx[4]*p[0];
		dx[5] = dx[5]*p[0];
	}

	bool ode_structure(DiffStructure &s) {
		/*s(0,x_indexode(2));
		s(0,x_indexode(3));
		s(0,p_indexode(0));

		s(1,x_indexode(2));
		s(1,x_indexode(3));
		s(1,p_indexode(0));

		s(2,u_indexode(0));
		s(2,p_indexode(0));

		s(3,x_indexode(2));
		s(3,x_indexode(4));
		s(3,p_indexode(0));

		s(4,u_indexode(1));
		s(4,p_indexode(0));

		s(5,u_indexode(0));
		s(5,u_indexode(1));*/
		return false;
	}

	bool ode_diff(DiffStructure& s, double t, const double* x, const double* u, const double* p) {
		/*s(0,x_indexode(2)) = cos(x[3]) * p[0];
		s(0,x_indexode(3)) = -x[2] * sin(x[3]) * p[0];

		s(1,x_indexode(2)) = sin(x[3]) * p[0];
		s(1,x_indexode(3)) = x[2] * cos(x[3]) * p[0];		

		s(2,u_indexode(0)) = p[0];

		s(3,x_indexode(2)) = sin(x[4]) * p[0] / 1; // b = 1
		s(3,x_indexode(4)) = x[2] * cos(x[4]) * p[0] / 1; // b = 1

		s(4,u_indexode(1)) = p[0];

		s(5,u_indexode(0)) = 2*u[0];
		s(5,u_indexode(1)) = 2*u[1];*/
		return false;
	}
	
	bool ode_diff_p(DiffStructure& s, double t, const double *x, const double *u, const double *p, int index) {
		/*s(0,p_indexode(0)) = x[2] * cos(x[3]);
		s(1,p_indexode(0)) = x[2] * sin(x[3]);
		s(2,p_indexode(0)) = u[0];
		s(3,p_indexode(0)) = x[2] * sin(x[4]) / 1; // b = 1
		s(4,p_indexode(0)) = u[1];*/
		return false;
	}
	
	void x_boundary(double *x_low, double *x_upp) {

		x_low[0] = 100; // min Geschw.
		x_upp[0] = 4000;

		x_low[1] = -M_PI; //Winkel
		x_upp[1] = M_PI;
		x_low[2] = -M_PI;
		x_upp[2] = M_PI;

		x_low[3] = 500; // Hoehe
		x_upp[3] = 100000;


		x_low[4] = -1.0; // lat
		x_upp[4] = 1.0;
		x_low[5] = -1.0; // lon
		x_upp[5] = 1.0;
		
	}
	
	void u_boundary(double *u_low, double *u_upp) {
		u_low[0] = 0.01; // Schranken aus PDF
		u_upp[0] = 0.18326;
		u_low[1] = -M_PI/2.0;
		u_upp[1] = M_PI/2.0;

	}

	void p_boundary(double *p_low, double *p_upp) {
		p_low[0] = 500; // sinnvolle Schranken fuer die Zeit
		p_upp[0] = 2000;
	}

	void var_boundary(double *x_low, double *x_upp) {

		// Anfangswerte
		x_low[x_index(0,0)] = x_upp[x_index(0,0)] = 2150.54529;
		x_low[x_index(0,1)] = x_upp[x_index(0,1)] = 0.152018177;
		x_low[x_index(0,2)] = x_upp[x_index(0,2)] = 2.2689279889;//95.0*M_PI/180.0;
		x_low[x_index(0,3)] = x_upp[x_index(0,3)] = 33900.0;
		x_low[x_index(0,4)] = x_upp[x_index(0,4)] = 0.8651416299;
		x_low[x_index(0,5)] = x_upp[x_index(0,5)] = 0.1980907302;
		
		// Endwert
		x_low[x_index(n_dis-1,3)] = x_upp[x_index(n_dis-1,3)] = 500.0;
	}

	void rand(double *r) {
		
	}
	
	void rand_boundary(double *r_low, double *r_upp) {
		
	}

};



/////////////////////////////////////////////////////////////////////////////

int main(int argv, char* argc[]) {

	TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv,argc);
	TWfolder folder(&twparameter,0);

	NotlandePhase ph(twparameter.NDIS);
	folder.Add(&ph);

	Viewer *viewer = 0;
	if (twparameter.PLOT) viewer = new Viewer(&twparameter);
	
	folder.Init();
	folder.Init(viewer);
	
	// Optional: Startschätzung
	/*
	for (int i=0;i<twparameter.NDIS;i++) {
		double z = i/(twparameter.NDIS-1.);
		ph.X[ph.x_index(i,0)] = 2000-z*1000;
		/*ph.X[ph.x_index(i,1)] = 2000;
		ph.X[ph.x_index(i,2)] = 2000;
		ph.X[ph.x_index(i,3)] = 2000;
		ph.X[ph.x_index(i,4)] = 2000;
		ph.X[ph.x_index(i,5)] = 2000;
* /
		ph.X[ph.u_index(i,0)] = 0.1;
		ph.X[ph.u_index(i,1)] = 0.8;
	}
	
	
	ph.Integrate(twparameter.butchertableau);
*/
	// Ende Optional
	
	folder.Loop(2);


	delete viewer;
	
	return 0;
}

