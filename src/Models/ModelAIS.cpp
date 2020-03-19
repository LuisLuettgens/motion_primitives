#include "ModelAIS.h"
using namespace std;



ModelAIS::ModelAIS ( string xmlparams ) :BaseModel ( xmlparams ), xmlparams(xmlparams) {
	
	
	nbrP = get_nbrP(  );
	nbrS = get_nbrS(  );
	nbrC = get_nbrC(  );
	
	get_paramAddr( paramAddr );
	loadParams (  );
}

void ModelAIS::ode ( double *dx, const double *x, const double *ctrl, const double *p, const double *q ) {
	x0        = x[0];
	y0        = x[1];
	psi       = x[2];
	u         = x[3];
	v         = x[4];
	r         = x[5];

	dx[0] = u * cos ( psi ) - v * sin ( psi );
	dx[1] = u * sin ( psi ) + v * cos ( psi );
	dx[2] = r;
	dx[3] = 0;
	dx[4] = 0;
	dx[5] = 0;
}

void ModelAIS::ode_structure ( tw::DiffStructure &s, int N ) {
	s ( N+0, N+2 );
	s ( N+0, N+3 );
	s ( N+0, N+4 );
	s ( N+1, N+2 );
	s ( N+1, N+3 );
	s ( N+1, N+4 );
	s ( N+2, N+5 );
}

void ModelAIS::loadParams (  ) {
	XMLParser p;
	XMLNode *n = p.Parse ( xmlparams )->GetFirstChild ( "PARAMS" );
	XMLNode *xml = n->GetFirstChild ( "DOUBLE" );
	
	vector<string> paramNames;
	get_paramNames( paramNames );
	
	while ( xml ) {
		string s ( xml->GetAttribute ( "name" ) );
		for( int i=0; i<nbrP; i++) {
			if ( s == paramNames[i] ) *paramAddr[i] = ToDouble ( xml->GetText() ); 
		}
		xml = n->GetNextChild ( "DOUBLE" );
	}
	init();
}

void ModelAIS::updateParams ( const double *p ) {
	XMLParser parser;
	XMLNode *n = parser.Parse ( xmlparams )->GetFirstChild ( "PARAMS" );
	XMLNode *xml = n->GetFirstChild ( "DOUBLE" );
	
	vector<string> paramNames;
	get_paramNames( paramNames );
	
	int j=0;
	while ( xml ) {
		string s1 ( xml->GetAttribute ( "name" ) );
		string s2 ( xml->GetAttribute ( "idfy" ) );
		
		for( int i=0; i<nbrP; i++) {
			if ( s1 == paramNames[i] and s2 == "TRUE" ) { *paramAddr[i] = p[j]; j++; };
		}
		xml = n->GetNextChild ( "DOUBLE" );
	}
	init();
}

void ModelAIS::init() {
}

void ModelAIS::get_bounds ( double *low, double *upp ) {
// 	string xmlparams = "ParamXML/s1v1_m1v3_p1.xml";
	XMLParser p;
	XMLNode *n1 = p.Parse ( xmlparams )->GetFirstChild ( "STATES" );
	XMLNode *xml1 = n1->GetFirstChild ( "DOUBLE" );
	
	int k=0;
	while ( xml1 ) {
		low[k] = ToDouble ( xml1->GetAttribute ( "low" ) );
		upp[k] = ToDouble ( xml1->GetAttribute ( "upp" ) );
		k++;
		xml1 = n1->GetNextChild ( "DOUBLE" );
	}
	
	XMLNode *n2 = p.Parse ( xmlparams )->GetFirstChild ( "CONTROLS" );
	XMLNode *xml2 = n2->GetFirstChild ( "DOUBLE" );
	
	while ( xml2 ) {
		low[k] = ToDouble ( xml2->GetAttribute ( "low" ) );
		upp[k] = ToDouble ( xml2->GetAttribute ( "upp" ) );
		k++;
		xml2 = n2->GetNextChild ( "DOUBLE" );
	}
}

void ModelAIS::get_param_bounds ( double *low, double *upp ) {
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

void ModelAIS::get_diff_bounds ( double *low, double *upp ) {
}


int ModelAIS::get_nbrS(  ) {
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

int ModelAIS::get_nbrC(  ) {
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

int ModelAIS::get_nbrP(  ) {
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

int ModelAIS::get_nbrIdParams(  ) {
	int nbrParams = 0;
	
	XMLParser p;
	XMLNode *n = p.Parse ( xmlparams )->GetFirstChild ( "PARAMS" );
	XMLNode *xml = n->GetFirstChild ( "DOUBLE" );
	
	while ( xml ) {
		string s ( xml->GetAttribute ( "idfy" ) );
		if ( s == "TRUE" ) nbrParams++;
		xml = n->GetNextChild ( "DOUBLE" );
	}
	return nbrParams;
}

void ModelAIS::get_paramNames( vector<string> &paramNames ) {
	XMLParser p;
	XMLNode *n = p.Parse ( xmlparams )->GetFirstChild ( "PARAMS" );
	XMLNode *xml = n->GetFirstChild ( "DOUBLE" );
	
	paramNames.resize( get_nbrP(  ) );
    int i=0;
	while ( xml ) {
		paramNames[i] = ( xml->GetAttribute ( "name" ) );
		i++;
		xml = n->GetNextChild ( "DOUBLE" );
	}
}

void ModelAIS::get_paramIdfy   ( vector<bool> &paramIdfy  ) {}
void ModelAIS::get_paramScales ( vector<double> &paramScales ) {}
void ModelAIS::get_initialS ( vector<double> &initialS ) {}

void ModelAIS::get_paramAddr( vector<double*> &paramAddr ) {

	paramAddr.resize( nbrP );
	
	vector<double*> pp = {
		&L,
		&BOA,
	};

	for(int i=0; i<nbrP; i++) paramAddr[i] = pp[i];
}












