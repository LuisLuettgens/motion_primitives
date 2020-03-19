#include "BaseModel.h"
using namespace std;

#ifndef MODEL1V2_H
#define MODEL1V2_H


class Model1v2: public BaseModel {
	public:

		double sign ( double x ) {
			return ( x >= 0 ) ? 1 : -1;
		}


		//Some constants
		double g; 				//Acceleration due to gravity (m/s^2)
		double rho; 				//denisty of water (kg/m^3)
		double rhoA;				//density of air (kg/m^3)

		// Ship main dimensions
		double L;				//Length between perpendiculars (m)
		double LWL;				//Length of waterline (m)
		double LOA;				//Length overall (m)
		double BOA;				//Beam (m)
		double dfor;				//Draft forward (m)
		double daft;				//Draft aft (m)
		double AF;				//Frontal area (m^2)
		double AL;				//Lateral area (m^2)
		double Cb;				//Block coefficient (-)
		double xG;				//Longitudinal center of buoyancy (LCB) forward of midship (m)
		double Rg;				//Radius of gyration (m)

		double l;				//Half-length (m)
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

		//Wind parameters
		double CX0;
		double CY0;
		double xA0;

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

		//Thruster paramters
		double ah;
		double bh;
		double xT;

		double aTY2;
		double aTY3;
		double bTY2;
		double bTY3;
		double aTN2;
		double aTN3;
		double bTN2;
		double bTN3;

		double kYm;
		double kYp;
		double kNm;
		double kNp;

		double uTYrel1;
		double uTYrel2;
		double uTYrel3;
		double uTYrel4;

		double uTNrel1;
		double uTNrel2;
		double uTNrel3;
		double uTNrel4;

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

		//Wind parameters
		double XA;
		double YA;
		double NA;

		double CXA;
		double CYA;
		double CNA;

		double AWx0;
		double AWy0;


		//Propellor1 effects
		double XP1;  		//Force  associated with "propeller" effects (x-coordinate) (N)
		double YP1;  		//Force  associated with "propeller" effects (y-coordinate) (N)
		double NP1;  		//Moment associated with "propeller" effects (about z-axis) (Nm)

		double cP1;		//Circumferential velocity of propeller blade at 0.7 radius (m/s)
		double uP1;		//Speed of advance of propeller (m/s)
		double eps1; 		//Advance angle (rad): tan(eps) = up/cp
		double CT1; 		//Propeller thrust coefficient (-)
		double T1; 		//Propeller thrust (N)
		double CQ1;
		double Q1;
		double QE1;
		double QEmax1;
		double YPT1;
		double NPT1;
		double RotS1;

		//Propellor2 effects
		double XP2;  		//Force  associated with "propeller" effects (x-coordinate) (N)
		double YP2;  		//Force  associated with "propeller" effects (y-coordinate) (N)
		double NP2;  		//Moment associated with "propeller" effects (about z-axis) (Nm)

		double cP2;		//Circumferential velocity of propeller blade at 0.7 radius (m/s)
		double uP2;		//Speed of advance of propeller (m/s)
		double eps2; 		//Advance angle (rad): tan(eps) = up/cp
		double CT2; 		//Propeller thrust coefficient (-)
		double T2; 		//Propeller thrust (N)
		double CQ2;
		double Q2;
		double QE2;
		double QEmax2;
		double YPT2;
		double NPT2;
		double RotS2;

		double xP;
		double yP;
		double YPTp;
		double YPTm;
		double NPTp;
		double NPTm;


		//Rudder1 effects
		double XR1;  		//Force  associated with "rudder" effects (x-coordinate) (N)
		double YR1;  		//Force  associated with "rudder" effects (y-coordinate) (N)
		double NR1;  		//Moment associated with "rudder" effects (about z-axis) (Nm)

		double uR1; 		//Longitudinal velocity at the rudder outside of the propeller slipstream (m/s)
		double vR1; 		//Transverse velocity at the rudder outside of the propeller slipstream (m/s)
		double uRP1; 		//Longitudinal velocity at the rudder in the propeller slipstream (m/s)
		double uRmean1; 	//Average longitudinal velocity at the rudder (m/s)
		double deltae1; 	//Effective angle of attack (rad)
		double betaR1; 	//Local drift angle (rad)
		double CLR01;
		double CDR01;
		double CLR1; 		//Rudder lift coefficient (-)
		double CDR1; 		//Rudder drag coefficient (-)

		//Rudder2 effects
		double XR2;  		//Force  associated with "rudder" effects (x-coordinate) (N)
		double YR2;  		//Force  associated with "rudder" effects (y-coordinate) (N)
		double NR2;  		//Moment associated with "rudder" effects (about z-axis) (Nm)

		double uR2; 		//Longitudinal velocity at the rudder outside of the propeller slipstream (m/s)
		double vR2; 		//Transverse velocity at the rudder outside of the propeller slipstream (m/s)
		double uRP2; 		//Longitudinal velocity at the rudder in the propeller slipstream (m/s)
		double uRmean2; 	//Average longitudinal velocity at the rudder (m/s)
		double deltae2; 	//Effective angle of attack (rad)
		double betaR2; 	//Local drift angle (rad)
		double CLR02;
		double CDR02;
		double CLR2; 		//Rudder lift coefficient (-)
		double CDR2; 		//Rudder drag coefficient (-)

		double xR; 		//Longitudinal position of rudder forward of midship (m)

		//Thruster
		double YT;
		double NT;
		double TT;

		double kY;
		double kN;

		double kh;

		double wjet;
		double uTrel;

		double aTY0;
		double aTY1;
		double bTY0;
		double bTY1;
		double aTN0;
		double aTN1;
		double bTN0;
		double bTN1;

		double TTmax;




		double m11, m22, m23, m32, m33;
		double det;
		double dv;
		double dr;


		Model1v2 ( string xmlparams );

		void ode ( double *dx, const double *x, const double *ctrl, const double *p, const double *q );

		void ode_structure ( tw::DiffStructure &s, int N );

		void loadParams ( string xmlparams );

		void init();

		void get_bounds ( double *low, double *upp );
		void get_diff_bounds ( double *low, double *upp );

		int get_nbrStates();
		int get_nbrControls();
};

#endif
