#include "Model1v4.h"
using namespace std;


	
double asign(double x) {
	return tanh(10000*x);
}

Model1v4::Model1v4 ( string xmlparams ) :BaseModel ( xmlparams ), xmlparams(xmlparams) {
	
	nbrP = get_nbrP(  );
	nbrS = get_nbrS(  );
	nbrC = get_nbrC(  );
	
	get_paramAddr( paramAddr );
	get_paramIdfy( paramIdfy );
	get_paramScales( paramScales );
	get_initialS( initialS );
	loadParams (  );
}

void Model1v4::ode ( double *dx, const double *x, const double *ctrl, const double *p, const double *q ) {
	l     = L / 2; 			//Half-length (m)
	d     = ( dfor + daft ) / 2; 		//Draft (m)
	AP    = M_PI * pow ( DP/2,2 );
	xR = ndxR * L;
	
	xP = ndxP * L;
	yP = ndyP * L;
	
	
	YPTp = ndYPTp;
	YPTm = ndYPTm;
	NPTp = ndNPTp * L;
	NPTm = ndNPTm * L;
	
	HR    = AR/ChordR;
// 	eta   = DP/HR;
	
// 	delta.resize(1);
// 	n.resize(1);
	
	x0        = x[0];
	y0        = x[1];
	psi       = x[2];
	u         = x[3];
	v         = x[4];
	r         = x[5];
// 	delta[0]  = x[6];
// 	n[0]      = x[7];
	
// 	deltaC.resize(1);
// 	EOT.resize(1);
// 	PP.resize(1);
	delta.resize(1);
	n.resize(1);
	PP.resize(1);

	X         = ctrl[0];
	Y         = ctrl[1];
	N         = ctrl[2];
	delta[0]  = ctrl[3];
	n[0]      = ctrl[4];
// 	deltaC[0] = ctrl[3];
// 	EOT[0]    = ctrl[4];
	PP[0]     = ctrl[5];
	
	DFT  = q[0];
	SET  = q[1];
	TWS  = q[2];
	TWA  = q[3];
	TWD  = q[4];
	AWS  = q[5];
	AWA  = q[6];
	AWD  = q[7];
	h    = q[8];
	
	updateParams( p );

	
	dx[0] = u * cos ( psi ) - v * sin ( psi );
	dx[1] = u * sin ( psi ) + v * cos ( psi );
	dx[2] = r;
	
// 	m = 2.272e+7;
	
	
	//Hull longitudinal resistance
	XHR = -m* ( XHRu/sqrt(L/g) * u + XHRuau/L * asign(u) * u*u + XHRuuu*sqrt(L/g)/(L*L) * u*u*u );
	
		//Hull lifting effects
	vs = k*r*l*asign ( u );
	sqroot = sqrt ( u*u + ( vs-v ) * ( vs-v ) );

	if ( sqroot < 1e-15 ) {
		XHL = 0;
		YHL = 0;
		NHL = 0;
	} else {
		XHL =  rho/2*L*  d * u* ( b1p*vs-b1*v ) /sqroot * ( ( vs-v ) - u*u* ( b2p*vs-b2*v ) / ( sqroot*sqroot ) );
		YHL =  rho/2*L*  d * u*u* ( b1p*vs-b1*v ) /sqroot * ( 1 + ( b2p*vs-b2*v ) * ( vs-v ) / ( sqroot*sqroot ) );
		NHL = -rho/2*L*L*d * asign(u)*u*u * ( b3p*vs-b3*v ) /sqroot * ( 1 + ( b2p*vs-b2*v ) * ( vs-v ) / ( sqroot*sqroot ) );
	}

	//Hull cross flow
	if ( v*v >= r*r * l*l ) {
		YHC = -rho/2.*L*  d * ( ( a0+a8/9. ) *pow ( v,2 ) +
								( a0/3.+a8/11. ) *pow ( r*l,2 ) +
								2.* ( a7/9. + a9/11. ) *v*r*l ) * asign ( v );

		NHC = -rho/2.*L*L*d * ( ( a7/18.+a9/22. ) *pow ( v,2 ) +
								( a7/22.+a9/26. ) *pow ( r*l,2 ) +
								( a0/3. + a8/11. ) *v*r*l ) * asign ( v );

	} else {
		YHC = -rho/2.*L*  d * ( a0* ( r*l + v*v/ ( 3.*r*l ) ) *v +
								a7/2.* ( v*v/4. + r*r*l*l/5. - v*v/180.*pow ( v/ ( r*l ),8 ) ) +
								a8/5.* ( r*l + v/99.*pow ( v/ ( r*l ),9 ) ) *v +
								a9/2.* ( v*v/5. + r*r*l*l/6. - v*v/330.*pow ( v/ ( r*l ),10 ) ) ) * asign ( r );

		NHC = -rho/2.*L*L*d * ( a0/4.* ( v*v + r*r*l*l/2. - pow ( v,4 ) / ( 6.*r*r*l*l ) ) +
								a7/10.* ( r*l + v/99.*pow ( v/ ( r*l ),9 ) ) *v +
								a8/4.* ( v*v/5. + r*r*l*l/6. - v*v/330.*pow ( v/ ( r*l ),10 ) ) +
								a9/12.* ( r*l + v/143.*pow ( v/ ( r*l ),11 ) ) *v ) * asign ( r );
	}
	
	//Propellor
	
	cP = 0.7 * M_PI * n[0] * DP * PP[0];
	uP = ( 1 - omP ) * u;

	if ( abs ( uP ) < 1e-15 and abs ( cP ) < 1e-15 ) {
		T = 0;
		Q = 0;
	} else {
		eps = atan2 ( uP,cP );

		if ( n[0] < 0 ) {
			YPT =  YPTm;
			NPT =  NPTm;
		} else {
			YPT = YPTp;
			NPT = NPTp;
		};

		CT = prT1*sin(2*(eps-eps_shift))*cos(eps) + prT2*sin(eps) + prT3*cos(2*(eps-eps_shift2))*sin(eps);
		T = rho/2.*AP* ( uP*uP + cP*cP ) *CT;
		
		
// 		CQ = 0;
// 		for ( int j=0; j<21; j++ ) CQ += AKQ[j]*cos ( j*(eps-eps_shift2) ) + BKQ[j]*sin ( j*(eps-eps_shift2) );
// 		Q = rho/2.*AP*DP* ( uP*uP + cP*cP ) *CQ;
	}

	XP1 = ( 1 - tP ) * T;
	YP1 = +RotS * YPT * T;
	NP1 = +RotS * NPT * T - yP*T;
	
	XP2 = ( 1 - tP ) * T;
	YP2 = -RotS * YPT * T;
	NP2 = -RotS * NPT * T  + yP*T;
	

// 	//Rudder
	T = Tfac*T;

	uR = ( 1 - omR ) * u;
	vR = ( v + r * xR ) * kHR;

	uRP = uP + asign ( T ) * ( ( kPR - 1./2 ) * asign ( T ) + 1./2 ) * ( sqrt ( pow ( uP,2 ) + 2*abs ( T ) / ( rho*AP ) ) - asign ( u ) *uP );

	if ( uRP*uP >= 0 ) {
		uRmean = asign ( uP ) *sqrt ( eta*uRP*uRP + ( 1-eta ) *uP*uP );
	} else {
		uRmean = asign ( uRP* ( eta*uRP*uRP - ( 1-eta ) *uP*uP ) ) *sqrt ( asign ( eta*uRP*uRP - ( 1-eta ) *uP*uP ) * ( eta*uRP*uRP - ( 1-eta ) *uP*uP ) );
	}
	
	betaR = atan2 ( -vR, uRmean );
	deltae = delta[0] + betaR;
	deltae = atan2 ( sin ( deltae ),cos ( deltae ) );

	CLR = li1*sin(2*deltae) + li2*sin(4*deltae);
	CDR = dr1*(1-cos(deltae));

	XR = rho / 2 * AR * ( pow ( uRmean,2 ) + pow ( vR,2 ) ) * ( CLR * sin ( betaR ) - CDR * cos ( betaR ) );
	YR = rho / 2 * AR * ( pow ( uRmean,2 ) + pow ( vR,2 ) ) * ( CLR * cos ( betaR ) + CDR * sin ( betaR ) );
	NR = YR * xR;

	//Wind resistance
	if ( std::isnan ( AWS ) ) {
		AWx0 = dx[0] + TWS*cos ( TWD );
		AWy0 = dx[1] + TWS*sin ( TWD );
		AWS = sqrt ( AWx0*AWx0 + AWy0*AWy0 );
		AWA = atan2 ( AWy0, AWx0 ) - psi;
	}
	CXA = -CX0*cos ( AWA );
	CYA = -CY0*sin ( AWA );
	CNA = -CY0*sin ( AWA ) * ( xA0/LOA + 1./4 - abs ( atan2 ( sin ( AWA ),cos ( AWA ) ) ) / ( 2*M_PI ) );
	XA =  rhoA/2*AWS*AWS*AF    *CXA;
	YA =  rhoA/2*AWS*AWS*AL    *CYA;
	NA =  rhoA/2*AWS*AWS*AL*LOA*CNA;
	
// 
// 	QE = c1*EOT[0]*EOT[0] + c2*n[0]*n[0] + c3*EOT[0]*n[0] + c4*EOT[0] + c5*n[0] + c6 + c7*EOT[0]*EOT[0]*EOT[0] + c8*n[0]*n[0]*n[0];

	double Xtotal = XHR+XHL+XP1+XP2+2*XR;//+XA;//X+
	double Ytotal = YHL+YHC+YP1+YP2+2*YR;//+YA;//Y+
	double Ntotal = NHL+NHC+NP1+NP2+2*NR;//+NA;//N+
	
	Xtotal =0.1*Xtotal;
	Ytotal =0.1*Ytotal;
	Ntotal =0.1*Ntotal;
	
	dx[3] =   (aa/m*Xtotal  + bb*v*r + L*cc*r*r + jj/L*v*v );
	dx[4] =   (dd/m*Ytotal  + ee/(m*L)*Ntotal   + ff*u*r);
	dx[5] =   (gg/m*Ytotal  + hh/(m*L)*Ntotal   + ii*u*r + 0*kk/L*u*v);
// 	dx[6] =   1./LagR * ( deltaC[0] - delta[0] );
// 	dx[7] =   1./ ( 2*M_PI*IEP ) * ( QE - Q );
	
}

void Model1v4::init() {
	l     = L / 2; 			//Half-length (m)
	d     = ( dfor + daft ) / 2; 		//Draft (m)
	AP    = M_PI * pow ( DP/2,2 );
	xR = ndxR * L;
	
	xP = ndxP * L;
	yP = ndyP * L;
	
	YPTp = ndYPTp;
	YPTm = ndYPTm;
	NPTp = ndNPTp * L;
	NPTm = ndNPTm * L;
	
	HR    = AR/ChordR;
	eta   = DP/HR;
}

void Model1v4::get_paramAddr( vector<double*> &paramAddr ) {

	paramAddr.resize( nbrP );
	
	vector<double*> pp = {
		&g,
		&rho,
		&rhoA,
		&L,
		&m,
		&LWL,
		&LOA,
		&AF,
		&AL,
		&BOA,
		&dfor,
		&daft,
		&Cb,
		&xG,
		&Rg,
		&aa,
		&bb,
		&cc,
		&dd,
		&ee,
		&ff,
		&gg,
		&hh,
		&ii,
		&jj,
		&kk,
		&XHRu,
		&XHRuau,
		&XHRuuu,
		&k,
		&b1,
		&b2,
		&b3,
		&b1p,
		&b2p,
		&b3p,
		&a0,
		&a7,
		&a8,
		&a9,
		&DP,
		&tP,
		&omP,
		&prT1,
		&prT2,
		&prT3,
		&RotS,
		&eps_shift,
		&eps_shift2,
		&ndxP,
		&ndyP,
		&ndYPTp,
		&ndYPTm,
		&ndNPTp,
		&ndNPTm,
		&AR,
		&ChordR,
		&omR,
		&li1,
		&li2,
		&dr1,
		&kHR,
		&kPR,
		&ndxR,
		&ndyR,
		&Tfac,
		&eta,
		&CX0,
		&CY0,
		&xA0,
		&LagR,
		&IEP,
		&c1,
		&c2,
		&c3,
		&c4,
		&c5,
		&c6,
		&c7,
		&c8,
		&c9,
		&c10,
	};
	for(int i=0; i<nbrP; i++) paramAddr[i] = pp[i];
}

void Model1v4::ode_structure ( tw::DiffStructure &s, int N ) {
	s ( N+0, N+2 );
	s ( N+0, N+3 );
	s ( N+0, N+4 );
	s ( N+1, N+2 );
	s ( N+1, N+3 );
	s ( N+1, N+4 );
	s ( N+2, N+5 );
	s ( N+3, N+3 );
	s ( N+3, N+4 );
	s ( N+3, N+5 ); 
	s ( N+4, N+3 );
	s ( N+4, N+4 );
	s ( N+4, N+5 ); 
	s ( N+5, N+3 );
	s ( N+5, N+4 );
	s ( N+5, N+5 ); 

	for ( int i=0; i<nbrC; i++ ) {
		s ( N+3, N+nbrS+i );
		s ( N+4, N+nbrS+i );
		s ( N+5, N+nbrS+i );
	}

}

void Model1v4::loadParams (  ) {
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

void Model1v4::updateParams ( const double *p ) {
	int j=0;
	for( int i=0; i<nbrP; i++) {
		if( paramIdfy[i] == true ) { *paramAddr[i] = paramScales[i] * p[j]; j++; };
	}
	init();
}

void Model1v4::get_bounds ( double *low, double *upp ) {
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

void Model1v4::get_param_bounds ( double *low, double *upp ) {
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

void Model1v4::get_diff_bounds ( double *low, double *upp ) {
}


int Model1v4::get_nbrIdParams(  ) {
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


void Model1v4::get_paramIdfy   ( vector<bool> &paramIdfy  ) {
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
	
void Model1v4::get_paramScales ( vector<double> &paramScales ) {
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

void Model1v4::get_initialS ( vector<double> &initialS ) {
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


int Model1v4::get_nbrP(  ) {
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

int Model1v4::get_nbrS(  ) {
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

int Model1v4::get_nbrC(  ) {
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







