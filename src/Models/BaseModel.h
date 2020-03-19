#ifndef BaseModel_H
#define BaseModel_H

#ifdef WIN32
#include "windows.h"
#endif
#include "TransWORHP.h"
#include "conversion.h"
#include "xmlio.h"
using namespace std;



class BaseModel {
	public:
// 		states
		double x0;
		double y0;
		double psi;
		double u;
		double v;
		double r;

// 		controls
		vector<double> deltaC;
		vector<double> delta;
		vector<double> EOT;
		vector<double> n;
		vector<double> PP;
		vector<double> kT;
		
		double X;
		double Y;
		double N;

		double DFT; //Current drift (m/s)
		double SET; //Current set (rad)

		double TWS; //True wind speed (m/s)
		double TWA; //True wind angle (rad)
		double TWD; //True wind direction (rad)

		double AWS; //Apparent wind speed (m/s)
		double AWA; //Apparent wind angle (rad)
		double AWD; //Apparent wind direction (rad)

		double h; //Depth (m)

		double L;
		double l;
		double BOA;
		
		int nbrP;
		int nbrS;
		int nbrC;
		vector<double*> paramAddr;
		vector<bool> paramIdfy;
		vector<double> paramScales;
		vector<double> initialS;

		BaseModel ( string xmlparams ) {}
		virtual void ode ( double *dx, const double *x,
						   const double *ctrl, const double *p, const double *q ) = 0;

		virtual void ode_structure ( tw::DiffStructure &s, int N ) = 0;

		virtual void loadParams   (  ) = 0;
		virtual void updateParams ( const double *p ) = 0;

		virtual void init() = 0;

		virtual void get_bounds ( double *low, double *upp ) = 0;
		virtual void get_diff_bounds ( double *low, double *upp ) = 0;
		virtual void get_param_bounds (  double *low, double *upp ) = 0;

		virtual int get_nbrS(  ) = 0;
		virtual int get_nbrC(  ) = 0;
		virtual int get_nbrP(  ) = 0;
		virtual void get_paramAddr( vector<double*> &paramAddr ) = 0;
		virtual void get_paramIdfy   ( vector<bool> &paramIdfy  ) = 0;
		virtual void get_paramScales ( vector<double> &paramScales ) = 0;
		virtual void get_initialS ( vector<double> &initialS ) = 0;
};

#endif
