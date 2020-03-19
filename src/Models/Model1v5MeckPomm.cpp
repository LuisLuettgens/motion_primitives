#include "Model1v5MeckPomm.h"
using namespace std;

double sin2diff(double up, double cp, double e, double del) {
	return 2*(up*cos(e)-cp*sin(e))*(cp*cos(e)+up*sin(e))/( (sqrt(up*up+cp*cp)+del)*(sqrt(up*up+cp*cp)+del) );
}

double cos2diff(double up, double cp, double e, double del) {
	return 1 - 2*(up*cos(e)-cp*sin(e))*(up*cos(e)-cp*sin(e))/( (sqrt(up*up+cp*cp)+del)*(sqrt(up*up+cp*cp)+del) );
}

double del = 0.0000001;


	
double Model1v5MeckPomm::asign(double x) {
	
	return tanh(10000*x);
}

Model1v5MeckPomm::Model1v5MeckPomm ( string xmlparams ) :BaseModel ( xmlparams ), xmlparams(xmlparams) {
	
	nbrP = get_nbrP(  );
	nbrS = get_nbrS(  );
	nbrC = get_nbrC(  );
	
	get_paramAddr  ( paramAddr   );
	get_paramIdfy  ( paramIdfy   );
	get_paramScales( paramScales );
	get_initialS( initialS );
	loadParams (  );
}

void Model1v5MeckPomm::ode ( double *dx, const double *x, const double *ctrl, const double *p, const double *q ) {
	
	Izz  = m*Rg*Rg*L*L;
	l    = L / 2; 
	AP   = M_PI * pow ( DP/2,2 );
	xR   = ndxR * L;
	xP   = ndxP * 28.95;
	yP   = ndyP * 28.95;
	YPTp = ndYPTp;
	YPTm = ndYPTm;
	NPTp = ndNPTp * L;
	NPTm = ndNPTm * L;
	HR   = AR/ChordR;

	x0   = x[0];
	y0   = x[1];
	psi  = x[2];
	u    = x[3];
	v    = x[4];
	r    = x[5];
	
	delta.resize(1);
	n    .resize(1);
	PP   .resize(1);

	delta[0]  = ctrl[0];
	n[0]      = ctrl[1];
	PP[0]     = 0.85;

// 	DFT  = q[0];
// 	SET  = q[1];
// 	TWS  = q[2];
// 	TWA  = q[3];
// 	TWD  = q[4];
// 	AWS  = q[5];
// 	AWA  = q[6];
// 	AWD  = q[7];
// 	h    = q[8];
	
	//Parameter ID only
    //updateParams( p );

	dx[0] = u * cos ( psi ) - v * sin ( psi );
	dx[1] = u * sin ( psi ) + v * cos ( psi );
	dx[2] = r;
	
	
	//Hull longitudinal resistance
	XHR = -m* ( XHRu/sqrt(L/g) * u + XHRuau/L * asign(u) * u*u + XHRuuu*sqrt(L/g)/(L*L) * u*u*u );
	
	//Hull lifting effects
	vs = k*r*l*asign ( u );
	sqroot = sqrt ( u*u + ( vs-v ) * ( vs-v ) ) + del;
	XHL =  rho/2*L*  d * u* ( b1p*vs-b1*v ) /sqroot * ( ( vs-v ) - u*u* ( b2p*vs-b2*v ) / ( sqroot*sqroot ) );
	YHL =  rho/2*L*  d * u*u* ( b1p*vs-b1*v ) /sqroot * ( 1 + ( b2p*vs-b2*v ) * ( vs-v ) / ( sqroot*sqroot ) );
	NHL = -rho/2*L*L*d * asign(u)*u*u * ( b3p*vs-b3*v ) /sqroot * ( 1 + ( b2p*vs-b2*v ) * ( vs-v ) / ( sqroot*sqroot ) );

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

	eps = atan2 ( uP,cP );
	
	YPT = (YPTp-YPTm)/2*asign(n[0])+(YPTp+YPTm)/2;
	NPT = (NPTp-NPTm)/2*asign(n[0])+(NPTp+NPTm)/2;

// 	CT = prT1*sin(2*(eps-eps_shift))*cos(eps) + prT2*sin(eps) + prT3*cos(2*(eps-eps_shift2))*sin(eps);
	CT = prT1*sin2diff(uP,cP,eps_shift,del)*cP/(sqrt(uP*uP+cP*cP)+del) + prT2*uP/(sqrt(uP*uP+cP*cP)+del) + prT3*cos2diff(uP,cP,eps_shift2,del)*uP/(sqrt(uP*uP+cP*cP)+del);
	T = rho/2.*AP* ( uP*uP + cP*cP ) *CT;
	Q=0;

	XP1 = ( 1 - tP ) * T;
	YP1 = +RotS * YPT * T;
	NP1 = +RotS * NPT * T - yP*T;


	//Rudder
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

// 	//Wind resistance
// 	if ( std::isnan ( AWS ) ) {
// 		AWx0 = dx[0] + TWS*cos ( TWD );
// 		AWy0 = dx[1] + TWS*sin ( TWD );
// 		AWS = sqrt ( AWx0*AWx0 + AWy0*AWy0 );
// 		AWA = atan2 ( AWy0, AWx0 ) - psi;
// 	}
// 	CXA = -CX0*cos ( AWA );
// 	CYA = -CY0*sin ( AWA );
// 	CNA = -CY0*sin ( AWA ) * ( xA0/LOA + 1./4 - abs ( atan2 ( sin ( AWA ),cos ( AWA ) ) ) / ( 2*M_PI ) );
// 	XA =  rhoA/2*AWS*AWS*AF    *CXA;
// 	YA =  rhoA/2*AWS*AWS*AL    *CYA;
// 	NA =  rhoA/2*AWS*AWS*AL*LOA*CNA;
	
	double Xtotal = XHR+XHL+2*XP1+2*XR/*+XA*/;	
	double Ytotal = YHL+YHC+YP1+2*YR/*+YA*/;
	double Ntotal = NHL+NHC+NP1+2*NR/*+NA*/;
	
	det1 = 1 - Xup;
	det2 = m*L*Nrp*Yvp - m*L*Nvp*Yrp + m*L*xG*Yrp - m*L*Nrp + m*xG*Nvp - m*xG*xG - Izz*Yvp + Izz;
	
	
	dx[3] =  1./det1*( 1./m*Xtotal + Xvv/L*v*v + (Xvr+1)*v*r + (L*Xrr+xG)*r*r );


	dx[4] =  1./det2*( (Izz-m*L*Nrp)/m*Ytotal + (L*Yrp-xG)*Ntotal - m*(L*Xup*Yrp-L*Yrp*Yvp-Xup*xG+Yvp*xG)/L*u*v + 
					(-m*L*Nrp*Xup + m*L*Nvp*Yrp - m*L*xG*Yrp+m*L*Nrp-m*xG*Nvp+m*xG*xG+Izz*Xup-Izz)*u*r );
	dx[5] =  1./det2*( (Nvp-xG)*Ytotal + (1-Yvp)*Ntotal + m*(Xup*Yvp-Yvp*Yvp-Xup+Yvp)/L*u*v + 
					m*(Nvp*Xup-Nvp*Yvp-xG*Xup+xG*Yvp)*u*r );
	
		
}

void Model1v5MeckPomm::init() {
}



void Model1v5MeckPomm::ode_structure ( tw::DiffStructure &s, int N ) {
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

void Model1v5MeckPomm::loadParams (  ) {
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

void Model1v5MeckPomm::updateParams ( const double *p ) {
	int j=0;
	for( int i=0; i<nbrP; i++) {
		if( paramIdfy[i] == true ) { *paramAddr[i] = paramScales[i] * p[j]; j++; };
	}
	init();
}

void Model1v5MeckPomm::get_bounds ( double *low, double *upp ) {
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

void Model1v5MeckPomm::get_diff_bounds ( double *low, double *upp ) {
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

void Model1v5MeckPomm::get_param_bounds ( double *low, double *upp ) {
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

void Model1v5MeckPomm::get_paramIdfy   ( vector<bool> &paramIdfy  ) {
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
	
void Model1v5MeckPomm::get_paramScales ( vector<double> &paramScales ) {
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

void Model1v5MeckPomm::get_initialS ( vector<double> &initialS ) {
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


int Model1v5MeckPomm::get_nbrP(  ) {
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

int Model1v5MeckPomm::get_nbrS(  ) {
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

int Model1v5MeckPomm::get_nbrC(  ) {
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

void Model1v5MeckPomm::get_paramAddr( vector<double*> &paramAddr ) {

	paramAddr.resize( nbrP );
	
	vector<double*> pp = {
		&g,
		&rho,
		&rhoA,
		&m,
		&xG,
		&Rg,
		&L,
		&d,
		&LOA,
		&BOA,
		&AF,
		&AL,
		&Xup,
		&Yvp,
		&Yrp,
		&Nvp,
		&Nrp,
		&Xvv,
		&Xvr,
		&Xrr,
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
		&Tfac,
		&eta,
		&CX0,
		&CY0,
		&xA0,
	};
	for(int i=0; i<nbrP; i++) paramAddr[i] = pp[i];
}





