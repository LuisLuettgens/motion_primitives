#include "BaseModel.h"
using namespace std;

#ifndef MODELAIS_H
#define MODELAIS_H


class ModelAIS: public BaseModel {
	public:

		double LOA;				//Length overall (m)
		
		string xmlparams;

		ModelAIS ( string xmlparams );

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
		int  get_nbrIdParams( );
		void get_paramNames ( vector<string > &paramNames );
		void get_paramAddr  ( vector<double*> &paramAddr  );
		void get_paramIdfy   ( vector<bool> &paramIdfy  );
		void get_paramScales ( vector<double> &paramScales );
		void get_initialS ( vector<double> &initialS );
};

#endif
