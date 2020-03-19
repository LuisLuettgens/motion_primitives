#include "BaseModel.h"
using namespace std;

#ifndef MODEL2V1_H
#define MODEL2V1_H


class Model2v1: public BaseModel {
	public:

		double sign ( double x ) {
			return ( x >= 0 ) ? 1 : -1;
		}


		//Some constants
		double g; 				//Acceleration due to gravity (m/s^2)
		double rho; 				//denisty of water (kg/m^3)

		// Ship main dimensions
		double L;				//Length between perpendiculars (m)
		double LWL;				//Length of waterline (m)
		double B;				//Beam (m)
		double dfor;				//Draft forward (m)
		double daft;				//Draft aft (m)
		double Cb;				//Block coefficient (-)
		double xG;				//Longitudinal center of buoyancy (LCB) forward of midship (m)
		double Rg;				//Radius of gyration (m)

		double d;				//Draft (m)
		double m;				//Mass (kg)
		double Izz;				//Moment of inertia



		//System parameters ####################################################################

		double ndXup;
		double ndXuu;
		double ndXvr;
		double ndXvv;
		double ndXcacdd;
		double ndXcacbd;
		double ndYvp;
		double ndYur;
		double ndYuv;
		double ndYvav;
		double ndYcacdad;
		double ndYcacbabad;
		double ndYT;
		double ndNrp;
		double ndNur;
		double ndNuv;
		double ndNanr;
		double ndNcacd;
		double ndNcacbabad;
		double ndNT;
		double ndTuu;
		double ndTun;
		double ndTnan;
		double ndCun;
		double ndCnn;

		double tt;


		//End of system parameters #############################################################


		double Xup;
		double Xuu;
		double Xvr;
		double Xvv;
		double Xcacdd;
		double Xcacbd;
		double Yvp;
		double Yur;
		double Yuv;
		double Yvav;
		double Ycacdad;
		double Ycacbabad;
		double YT;
		double Nrp;
		double Nur;
		double Nuv;
		double Nanr;
		double Ncacd;
		double Ncacbabad;
		double NT;
		double Tuu;
		double Tun;
		double Tnan;
		double Cun;
		double Cnn;


		double beta;
		double T;
		double C;


		Model2v1 ( string xmlparams );

		void ode ( double *dx, const double *x, const double *ctrl, const double *p, const double *q );

		void loadParams ( string xmlparams );

		void init();
};

#endif
