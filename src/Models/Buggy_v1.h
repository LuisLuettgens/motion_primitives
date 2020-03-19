#ifndef BUGGY_V1_H
#define BUGGY_V1_H

#include "BaseModel.h"


class Buggy_v1: public BaseModel {
	public:

		double m;
		double g;
		double RW;
		double K_M;
		double P_DC;
		double P_DV;
		double P_T;
		double P_V;
		double totalMass;
		double thetaPowertrain;
		double thetaTire;
		
		string xmlparams;

		Buggy_v1 ( string xmlparams );

		void ode ( double *dx, const double *x, const double *ctrl, const double *p, const double *q );

		void ode_structure ( tw::DiffStructure &s, int N );

		void loadParams   (  );
		void updateParams ( const double *p );

		void init();

		void get_bounds       ( double *low, double *upp    );
		void get_diff_bounds  ( double *low, double *upp    );
		void get_param_bounds ( double *low, double *upp    );

		int  get_nbrS         (  );
		int  get_nbrC         (  );
		int  get_nbrP         (  );
		void get_paramAddr    ( vector<double*> &paramAddr  );
		void get_paramIdfy    ( vector<bool> &paramIdfy     );
		void get_paramScales  ( vector<double> &paramScales );
		void get_initialS     ( vector<double> &initialS    );
};

#endif
