#include "BaseModel.h"
using namespace std;

#ifndef MODEL1V5_H
#define MODEL1V5_H


class Model1v5: public BaseModel {
	public:

		double asign(double x);
		double sign ( double x ) {
			return ( x >= 0 ) ? 1 : -1;
		}


		//Some constants
		double g; 				//Acceleration due to gravity (m/s^2)
		double rho; 				//denisty of water (kg/m^3)
		double rhoA;				//density of air (kg/m^3)

		// Ship main dimensions
//   double L;					//Length between perpendiculars (m)
		double LWL;				//Length of waterline (m)
		double LOA;				//Length overall (m)
//   double BOA;				//Beam (m)
		double dfor;				//Draft forward (m)
		double daft;				//Draft aft (m)
		double AF;				//Frontal area (m^2)
		double AL;				//Lateral area (m^2)
		double Cb;				//Block coefficient (-)
		double xG;				//Longitudinal center of buoyancy (LCB) forward of midship (m)
		double Rg;				//Radius of gyration (m)

//   double l;				//Half-length (m)
		double d;				//Draft (m)
		double m;				//Mass (kg)
		double Izz;				//Moment of inertia about z-axis (kg*m^2)

		double DP;				//Diameter of propellors (m)
		double AP;				//Area of propellors (m^2)

		double AR;				//Area of rudders (m^2)
		double ChordR;			//Chord length of rudders (m)
		double HR;				//Height of rudders (m)

		double hT;				//Depth of thruster's axis (m)
		double RT;				//Radius of thruster (m)
		double AT;				//Area of thruster (m^2)

		double eta;				//Ratio between "diameter of propeller" and "height of rudder" (-)


		//System parameters ####################################################################

		//Ideal fluid parameters (non-dimensional)
		double ndXup;
		double ndXrr;
		double ndYvp;
		double ndNrp;
		double ndXvr;
		double ndXvv;
		double ndYrp;
		double ndNvp;

		//Hull lifting parameters
		double b1;
		double b2;
		double b3;
		double b1p;
		double b2p;
		double b3p;
		double k;

		//Hull nonlinear resistance parameters
		double ndXHRu;
		double ndXHRuau;
		double ndXHRuuu;
		double a0;
		double a7;
		double a8;
		double a9;

		//Interaction factors
		double tP;
		double omP;
		double omR;

		//Propellor parameters
		double ndxP;
		double ndyP;
		double IEP;
		double prT[11];
		double prtT;
		double AKT[21];
		double BKT[21];
		double prQ[13];
		double prtQ;
		double AKQ[21];
		double BKQ[21];
		double ndYPTp;
		double ndYPTm;
		double ndNPTp;
		double ndNPTm;

		//Rudder parameters
		double LagR;
		double kHR;
		double kPR;
		double ndxR;
		double ndyR;
		double kLR;
		double kDR;
		double kNR;
		double kCLR;
		double kCDR;

		double lift0, lift1, lift2, lift3, lift4, lift5;
		double drag0, drag1, drag2, drag3, drag4, drag5;
		double lia0, lia1, lia2, lia3;
		double lib0, lib1, lib2, lib3;
		double dra0, dra1, dra2, dra3;
		double drb0, drb1, drb2, drb3;
		double drc0, drc1, drc2, drc3;
		double de0, de1, de2, de3, de4, de5;
		
		//Wind parameters
		double CX0;
		double CY0;
		double xA0;
		
		double XA;
		double YA;
		double NA;

		double CXA;
		double CYA;
		double CNA;

		double AWx0;
		double AWy0;
		
		double Tfac;
		double ppp1;
		double ppp2;
		double ppp3;
		double pEOT;
		double c1,c2,c3,c4,c5,c6,c7,c8,c9,c10;
		double cc1,cc2,cc3;
		double ccc1,ccc2,ccc3,ccc4,ccc5;
		double ddd1;
		
		double prT1, prT2, prT3;
		double prQ1, prQ2, prQ3;
		
		double li1,li2;
		double dr1;

		//End of system parameters #############################################################


		//Ideal fluid effects
		double XI; 		//Force  associated with "ideal fluid" effects (x-coordinate) (N)
		double YI; 		//Force  associated with "ideal fluid" effects (y-coordinate) (N)
		double NI; 		//Moment associated with "ideal fluid" effects (about z-axis) (Nm)

		double Xup;		//Ideal fluid parameters
		double Xrr;		//-
		double Yvp;		//-
		double Nrp;		//-
		double Xvr;		//-
		double Xvv;		//-
		double Yrp;		//-
		double Nvp;		//-


		//Hull lifting effects
		double XHL; 		//Force  associated with "hull lifting" effects (x-coordinate) (N)
		double YHL; 		//Force  associated with "hull lifting" effects (y-coordinate) (N)
		double NHL; 		//Moment associated with "hull lifting" effects (about z-axis) (Nm)

		double vs;
		double sqroot;

		//Hull longitudinal resistance
		double XHR;
		//Hull cross flow
		double YHC;
		double NHC;

		double XHRu;		//Hull resistance parameters
		double XHRuau;	//-
		double XHRuuu;	//-


		//Propellor effects
		double XP;  		//Force  associated with "propeller" effects (x-coordinate) (N)
		double YP;  		//Force  associated with "propeller" effects (y-coordinate) (N)
		double NP;  		//Moment associated with "propeller" effects (about z-axis) (Nm)
		
		double XP1;  		//Force  associated with "propeller" effects (x-coordinate) (N)
		double YP1;  		//Force  associated with "propeller" effects (y-coordinate) (N)
		double NP1;  		//Moment associated with "propeller" effects (about z-axis) (Nm)
		
		double XP2;  		//Force  associated with "propeller" effects (x-coordinate) (N)
		double YP2;  		//Force  associated with "propeller" effects (y-coordinate) (N)
		double NP2;  		//Moment associated with "propeller" effects (about z-axis) (Nm)

		double cP;		//Circumferential velocity of propeller blade at 0.7 radius (m/s)
		double uP;		//Speed of advance of propeller (m/s)
		double eps; 		//Advance angle (rad): tan(eps) = up/cp
		double eps_shift;
		double eps_shift2;
		double CT; 		//Propeller thrust coefficient (-)
		double T; 		//Propeller thrust (N)
		double CQ;
		double Q;
		double QE;
		double QEmax;
		double YPT;
		double NPT;
		double RotS;

		double xP;
		double yP;
		double YPTp;
		double YPTm;
		double NPTp;
		double NPTm;


		//Rudder effects
		double XR;  		//Force  associated with "rudder" effects (x-coordinate) (N)
		double YR;  		//Force  associated with "rudder" effects (y-coordinate) (N)
		double NR;  		//Moment associated with "rudder" effects (about z-axis) (Nm)

		double uR; 		//Longitudinal velocity at the rudder outside of the propeller slipstream (m/s)
		double vR; 		//Transverse velocity at the rudder outside of the propeller slipstream (m/s)
		double uRP; 		//Longitudinal velocity at the rudder in the propeller slipstream (m/s)
		double uRmean; 	//Average longitudinal velocity at the rudder (m/s)
		double deltae; 	//Effective angle of attack (rad)
		double betaR; 	//Local drift angle (rad)
		double CLR0;
		double CDR0;
		double CLR; 		//Rudder lift coefficient (-)
		double CDR; 		//Rudder drag coefficient (-)

		double xR; 		//Longitudinal position of rudder forward of midship (m)

        double aa;
		double bb;
		double cc;
		double dd;
		double ee;
		double ff;
		double gg;
		double hh;
		double ii;
		double jj;
		double kk;
		double ll, mm, nn, oo;
		

		double m11, m22, m23, m32, m33;
		double det;
		double dv;
		double dr;
		
		double det1, det2;
		
		string xmlparams;
		vector<string> paramNames;

		Model1v5 ( string xmlparams );

		void ode ( double *dx, const double *x, const double *ctrl, const double *p, const double *q );

		void ode_structure ( tw::DiffStructure &s, int N );

		void loadParams   (  );
		void updateParams ( const double *p );

		void init();

		void get_bounds ( double *low, double *upp );
		void get_diff_bounds ( double *low, double *upp );
		void get_param_bounds ( double *low, double *upp );

		int  get_nbrS       (  );
		int  get_nbrC       (  );
		int  get_nbrP       (  );
		int  get_nbrIdParams(  );
		void get_paramAddr   ( vector<double*> &paramAddr  );
		void get_paramIdfy   ( vector<bool> &paramIdfy  );
		void get_paramScales ( vector<double> &paramScales );
		void get_initialS ( vector<double> &initialS );
};

#endif
