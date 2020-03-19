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

//double phasenPlot (int &len, int &ndgl, int &nsteuer, double *t,double *x, double *u, int &i, int &index);

class NotlandePhase : public TransWorhp {
public:

	double TMAX;
	double ISP;
	double FF;
	double CD0;
	double K;
	double RHO0;
	double G0;
	double R0;
	double HSKAL;
	double BETA;
	double W;
	

	NotlandePhase(int dis) : TransWorhp("Notlande", dis,6,2,0,0,0) {
		freetime = 0;

	// KONSTANTEN
	TMAX  = 1.2e+6;
	ISP   = 472.0;
	FF    = 305.0;         // Alternativ: 135.0D0
	CD0   = 0.017;
	K     = 2.0;
	RHO0  = 1.249512 + 1.249512*1/100.0; //1.249512 + 1.249512*p(1)/100.0;
	G0    = 9.80665;
	R0    = 6.371e+6;
	HSKAL = 6900.0;
	BETA  = 1.0/HSKAL;
	W     = 7.27e-5;

	}


	void OpenWindows(Viewer *viewer) {
		
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
	
	}

	
	double obj() {

		double c1 = 1;
		double c2 = 1;		

		return -c1*(x(n_dis-1,4)*x(0, 4))*(x(n_dis-1,4)*x(0, 4)) - c2*(x(n_dis-1,5)*x(0, 5))*(x(n_dis-1,5)*x(0, 5));
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
		double R     = R0+x[3];
		double GH    = G0*R0*R0/(R*R);
		double RHOH  = RHO0*exp(-BETA*x[3]);
		double QH    = 0.5*RHOH*x[0]*x[0];
		double CD    = CD0+K*u[0]*u[0];

		double HILF1 = QH*FF;
		double D     = HILF1*CD;
		double L     = HILF1*u[0];

		// ABKUERZUNGEN
		double COSU3   = cos(0);
		double SINU3   = sin(0);
		double COSX2   = cos(x[1]);
		double SINX2   = sin(x[1]);
		double COSX5   = cos(x[4]);
		double SINX5   = sin(x[4]);
		double COSX3   = cos(x[2]);
		double SINX3   = sin(x[2]);    
		double WT2     = 2.0*W;
		double W2COSX5 = W*W*COSX5;
		
 		double X7 = 115000; // massse
		
		dx[0] = -GH*SINX2 + R*W2COSX5*(SINX2*COSX5-COSX2*SINX3*SINX5);
		dx[1] = 0.0*cos(u[1])/(X7*x[0]) - (GH/x[0]-x[0]/R)*COSX2 + WT2*COSX3*COSX5 + W2COSX5*(SINX2*SINX3*SINX5+COSX2*COSX5)*R/x[0];
		dx[2] = 0.0*sin(u[1])/(X7*x[0]*COSX2) - COSX2*COSX3*tan(x[4])*x[0]/R + WT2*(SINX3*COSX5*tan(x[1])-SINX5) - W2COSX5*SINX5*COSX3*R/(x[0]*COSX2);
		dx[3] = x[0]*SINX2,
		dx[4] = COSX2*SINX3*x[0]/R;
		dx[5] = COSX2*COSX3*x[0]/(R*COSX5);


		//cout << (X7*x[0]) << endl;
		
		//cout << dx[0] << " " << dx[1] << " " << dx[2] << " " << dx[3] << " " << dx[4] << " " << dx[5] << endl;
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

		/*
		 BL(1) = -1.0d0
      BU(1) = 1.0d0
      BL(2) = -1.0d0
      BU(2) = 1.0d0
		*/

		x_low[0] = 100.0;
		x_low[1] = -1.0;
		x_upp[1] = 1.0;
	}
	
	void u_boundary(double *u_low, double *u_upp) {
		/*u_low[0] = -10;
		u_upp[0] = 10;
		u_low[1] = -M_PI/2; // sigma_max
		u_upp[1] = M_PI/2;*/
	}

	void p_boundary(double *p_low, double *p_upp) {
		
	}

	void var_boundary(double *x_low, double *x_upp) {

/*
		 AWX(1) = 2150.54529D0             ! 6.8 Mach Geschwindigkeit
      AWX(2) = 0.152018177d0            ! 8.71 Grad Bahnneigungswinkel
c      AWX(3) = 1.570796327d0            ! 90 Grad Kurswinkel
C      AWX(3) = 130.0d0*3.1415926d0/180.0d0
       AWX(3) = 95.0d0*3.1415926d0/180.0d0
       AWX(4) = 33900.0d0   +p(3)         ! 33.9 km Hoehe
c      AWX(5) = 0.287979326d0            ! 16.5 Grad geographische Breite
c      AWX(5) = 3.141592654d0            
c      AWX(6) = 0.060039326d0            ! 3.44 Grad geographische Laenge
      AWX(5) = DLAT          ! Breitengrad Bayreuth  
      AWX(6) = DLON          ! Laengengrad Bayreuth
      AWX(7) = 115000.0d0               ! Kg Masse
c      AWX(7) = 0.0d0
      AWX(8) = unknown(1)
*/

		x_low[x_index(0,0)] = x_upp[x_index(0,0)] = 2150.54529;
		x_low[x_index(0,1)] = x_upp[x_index(0,1)] = 0.152018177;
		x_low[x_index(0,2)] = x_upp[x_index(0,2)] = 95.0*3.1415926/180.0;
		x_low[x_index(0,3)] = x_upp[x_index(0,3)] = 33900.0;
		x_low[x_index(0,4)] = x_upp[x_index(0,4)] = 0.8651416299;
		x_low[x_index(0,5)] = x_upp[x_index(0,5)] = 0.1980907302;
		
		x_low[x_index(n_dis-1,3)] = x_upp[x_index(n_dis-1,3)] = 10.0;
		x_low[x_index(n_dis-1,4)] = x_upp[x_index(n_dis-1,4)] = 0.92692012505013;
		x_low[x_index(n_dis-1,5)] = x_upp[x_index(n_dis-1,5)] = 0.15452466458004;
	}

	void rand(double *r) {
		//r[0] = x(3,n_dis-1)-10.0;// + p(2) //! 500
		//r[0] = r[0]*1.0e-3;
	}
	
	void rand_boundary(double *r_low, double *r_upp) {
		r_upp[0] = 0;
		r_upp[1] = 0;
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
	folder.Loop();

	delete viewer;
	
	return 0;
}

