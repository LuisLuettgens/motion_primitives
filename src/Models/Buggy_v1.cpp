#include "Buggy_v1.h"
using namespace std;


double W(double v) {
	return (1.0 / (exp(-50.0*v) + 1.0) - 0.5)*2.0;
}


Buggy_v1::Buggy_v1 ( string xmlparams ) :BaseModel ( xmlparams ), xmlparams(xmlparams) {
	
	nbrP = get_nbrP(  );
	nbrS = get_nbrS(  );
	nbrC = get_nbrC(  );
	
	get_paramAddr  ( paramAddr   );
	get_paramIdfy  ( paramIdfy   );
	get_paramScales( paramScales );
	get_initialS( initialS );
	loadParams (  );
}

void Buggy_v1::ode ( double *dx, const double *x, const double *ctrl, const double *p, const double *q ) {
	
	m = totalMass + (2.0*thetaPowertrain + 4.0*thetaTire) / (RW*RW);

	x0   = x[0];
	y0   = x[1];
	psi  = x[2];
	u    = x[3];
	v    = x[4];
	r    = x[5];
	
	delta.resize(1);
	n    .resize(1);

	delta[0]  = ctrl[0];
	n[0]      = ctrl[1];
	
	dx[0] = x[3] * cos(x[2]);
	dx[1] = x[3] * sin(x[2]);
	dx[2] = (x[3] / L)*tan(delta[0]);
	dx[3] = ((n[0] * K_M*4.8 / RW) - (x[3] * P_V + W(x[3])*P_T*g*totalMass + x[3] * P_DV / (RW*RW) + W(x[3])*P_DC / RW)) / m;
	dx[4] = 0;
	dx[5] = 0;
}

void Buggy_v1::init() {
}



void Buggy_v1::ode_structure ( tw::DiffStructure &s, int N ) {
	s ( N+0, N+2 );
	s ( N+0, N+3 );
	s ( N+1, N+2 );
	s ( N+1, N+3 );
	s ( N+2, N+3 );
	s ( N+3, N+3 );
	
	s ( N+2, N+nbrS+0 );
	s ( N+3, N+nbrS+1 );
}

void Buggy_v1::loadParams (  ) {
	XMLParser p;
	XMLNode *n = p.Parse ( xmlparams )->GetFirstChild ( "PARAMS" );
	XMLNode *xml = n->GetFirstChild ( "DOUBLE" );
	
	int i=0;
	while ( xml ) {
		string s ( xml->GetAttribute ( "name" ) );
		string scale ( xml->GetAttribute ( "scale" ) );
		*paramAddr[i] = ToDouble( scale ) * ToDouble ( xml->GetText() ); 
		xml = n->GetNextChild ( "DOUBLE" );
		i++;
	}
	init();
}

void Buggy_v1::updateParams ( const double *p ) {
	int j=0;
	for( int i=0; i<nbrP; i++) {
		if( paramIdfy[i] == true ) { *paramAddr[i] = paramScales[i] * p[j]; j++; };
	}
	init();
}

void Buggy_v1::get_bounds ( double *low, double *upp ) {
	XMLParser p;
	
	int k=0;
	
	XMLNode *n1 = p.Parse ( xmlparams )->GetFirstChild ( "STATES" );
	XMLNode *xml1 = n1->GetFirstChild ( "DOUBLE" );
	while ( xml1 ) {
		low[k] = ToDouble ( xml1->GetAttribute ( "low" ) );
		upp[k] = ToDouble ( xml1->GetAttribute ( "upp" ) );
		xml1 = n1->GetNextChild ( "DOUBLE" );
		k++;
	}
	
	XMLNode *n2 = p.Parse ( xmlparams )->GetFirstChild ( "CONTROLS" );
	XMLNode *xml2 = n2->GetFirstChild ( "DOUBLE" );
	while ( xml2 ) {
		low[k] = ToDouble ( xml2->GetAttribute ( "low" ) );
		upp[k] = ToDouble ( xml2->GetAttribute ( "upp" ) );
		xml2 = n2->GetNextChild ( "DOUBLE" );
		k++;
	}
}

void Buggy_v1::get_diff_bounds ( double *low, double *upp ) {
	XMLParser p;
	
	int k=0;
	
	XMLNode *n2 = p.Parse ( xmlparams )->GetFirstChild ( "CONTROLS" );
	XMLNode *xml2 = n2->GetFirstChild ( "DOUBLE" );
	while ( xml2 ) {
		low[k] = ToDouble ( xml2->GetAttribute ( "dlow" ) );
		upp[k] = ToDouble ( xml2->GetAttribute ( "dupp" ) );
		xml2 = n2->GetNextChild ( "DOUBLE" );
		k++;
	}
}

void Buggy_v1::get_param_bounds ( double *low, double *upp ) {
	XMLParser p;
	XMLNode *n = p.Parse ( xmlparams )->GetFirstChild ( "PARAMS" );
	XMLNode *xml = n->GetFirstChild ( "DOUBLE" );
	
	int k=0;
	while ( xml ) {
		string s ( xml->GetAttribute ( "idfy" ) );
		if ( s == "TRUE" ) {
			low[k] = ToDouble ( xml->GetAttribute ( "low" ) );
			upp[k] = ToDouble ( xml->GetAttribute ( "upp" ) );
			k++;
		}
		xml = n->GetNextChild ( "DOUBLE" );
	}
}

void Buggy_v1::get_paramIdfy   ( vector<bool> &paramIdfy  ) {
	XMLParser p;
	XMLNode *n = p.Parse ( xmlparams )->GetFirstChild ( "PARAMS" );
	XMLNode *xml = n->GetFirstChild ( "DOUBLE" );
	
	paramIdfy.resize( get_nbrP(  ) );
    int i=0;
	while ( xml ) {
		paramIdfy[i] = xml->GetAttribute ( "idfy" ) == "TRUE" ? true : false;
		i++;
		xml = n->GetNextChild ( "DOUBLE" );
	}
}
	
void Buggy_v1::get_paramScales ( vector<double> &paramScales ) {
	XMLParser p;
	XMLNode *n = p.Parse ( xmlparams )->GetFirstChild ( "PARAMS" );
	XMLNode *xml = n->GetFirstChild ( "DOUBLE" );
	
	paramScales.resize( get_nbrP(  ) );
    int i=0;
	while ( xml ) {
		paramScales[i] = ToDouble( xml->GetAttribute ( "scale" ) );
		i++;
		xml = n->GetNextChild ( "DOUBLE" );
	}
}

void Buggy_v1::get_initialS ( vector<double> &initialS ) {
	XMLParser p;
	XMLNode *n = p.Parse ( xmlparams )->GetFirstChild ( "STATES" );
	XMLNode *xml = n->GetFirstChild ( "DOUBLE" );
	
	initialS.resize( get_nbrS(  ) );
    int i=0;
	while ( xml ) {
		initialS[i] = ToDouble( xml->GetAttribute ( "init" ) );
		i++;
		xml = n->GetNextChild ( "DOUBLE" );
	}
}


int Buggy_v1::get_nbrP(  ) {
	int nbrParams = 0;
	
	XMLParser p;
	XMLNode *n = p.Parse ( xmlparams )->GetFirstChild ( "PARAMS" );
	XMLNode *xml = n->GetFirstChild ( "DOUBLE" );
	
	while ( xml ) {
		nbrParams++;
		xml = n->GetNextChild ( "DOUBLE" );
	}
	return nbrParams;
}

int Buggy_v1::get_nbrS(  ) {
	int nbrStates = 0;
	
	XMLParser p;
	XMLNode *n = p.Parse ( xmlparams )->GetFirstChild ( "STATES" );
	XMLNode *xml = n->GetFirstChild ( "DOUBLE" );
	
	while ( xml ) {
		nbrStates++;
		xml = n->GetNextChild ( "DOUBLE" );
	}
	return nbrStates;
}

int Buggy_v1::get_nbrC(  ) {
	int nbrControls = 0;
	
	XMLParser p;
	XMLNode *n = p.Parse ( xmlparams )->GetFirstChild ( "CONTROLS" );
	XMLNode *xml = n->GetFirstChild ( "DOUBLE" );
	
	while ( xml ) {
		nbrControls++;
		xml = n->GetNextChild ( "DOUBLE" );
	}
	return nbrControls;
}

void Buggy_v1::get_paramAddr( vector<double*> &paramAddr ) {

	paramAddr.resize( nbrP );
	
	vector<double*> pp = {
		&g,
		&L,
		&RW,
		&K_M,
		&P_DC,
		&P_DV,
		&P_T,
		&P_V,
		&totalMass,
		&thetaPowertrain,
		&thetaTire,
	};
	for(int i=0; i<nbrP; i++) paramAddr[i] = pp[i];
}





