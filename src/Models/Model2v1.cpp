#include "Model2v1.h"
using namespace std;


Model2v1::Model2v1 ( string xmlparams ) :BaseModel ( xmlparams ) {

	loadParams ( xmlparams );
	init();

	delta.resize ( 1 );
	n.resize ( 1 );
}

void Model2v1::ode ( double *dx, const double *x, const double *ctrl, const double *p, const double *q ) {
	xo     = x[0];
	yo     = x[1];
	psi    = x[2];
	u      = x[3];
	v      = x[4];
	r      = x[5];

	delta[0]  = ctrl[0];
	n[0]      = ctrl[1];

	U_c     = q[0];
	psi_c   = q[1];
	U_tw     = q[2];
	psi_tw   = q[3];
	h      = q[4];


	beta = atan2 ( -v,u );

	T = Tuu*u*u + Tun*u*n[0] + Tnan*n[0]*abs ( n[0] );
	C = ( n[0] >= 0 ) ? Cun*u*n[0] + Cnn*n[0]*n[0] : 0;

	//Dynamics
	dx[0] = u * cos ( psi ) - v * sin ( psi ) + U_c*cos ( psi_c );
	dx[1] = u * sin ( psi ) + v * cos ( psi ) + U_c*sin ( psi_c );
	dx[2] = r;

	dx[3] = 1./ ( m-Xup )   * ( Xuu*u*u + ( m+Xvr ) *v*r + Xvv*v*v + Xcacdd*C*abs ( C ) *delta[0]*delta[0] +
								Xcacbd*C*abs ( C ) *beta*delta[0] + ( 1-tt ) *T );
	dx[4] = 1./ ( m-Yvp )   * ( ( Yur-m ) *u*r + Yuv*u*v + Yvav*v*abs ( v ) + Ycacdad*C*abs ( C ) *delta[0]*abs ( delta[0] ) +
								Ycacbabad*C*abs ( C ) *beta*abs ( beta ) *abs ( delta[0] ) + YT*T );
	dx[5] = 1./ ( Izz-Nrp ) * ( ( Nur-m*xG ) *u*r + Nuv*u*v + Nanr*abs ( v ) *r +
								Ncacd*C*abs ( C ) *delta[0] + Ncacbabad*C*abs ( C ) *beta*abs ( beta ) *abs ( delta[0] ) +
								NT*T );

}

void Model2v1::loadParams ( string xmlparams ) {
	XMLParser p;
	XMLNode *n = p.Parse ( xmlparams )->GetFirstChild ( "PARAMS" );
	XMLNode *xml = n->GetFirstChild ( "DOUBLE" );

	while ( xml ) {
		string s ( xml->GetAttribute ( "name" ) );

		if ( s == "L" )  L           = ToDouble ( xml->GetText() );

		if ( s == "LWL" )  LWL         = ToDouble ( xml->GetText() );

		if ( s == "B" )  B           = ToDouble ( xml->GetText() );

		if ( s == "dfor" )  dfor        = ToDouble ( xml->GetText() );

		if ( s == "daft" )  daft        = ToDouble ( xml->GetText() );

		if ( s == "Cb" )  Cb          = ToDouble ( xml->GetText() );

		if ( s == "xG" )  xG          = ToDouble ( xml->GetText() );

		if ( s == "Rg" )  Rg          = ToDouble ( xml->GetText() );

		if ( s == "ndXup" )  ndXup       = ToDouble ( xml->GetText() );

		if ( s == "ndXuu" )  ndXuu       = ToDouble ( xml->GetText() );

		if ( s == "ndXvr" )  ndXvr       = ToDouble ( xml->GetText() );

		if ( s == "ndXvv" )  ndXvv       = ToDouble ( xml->GetText() );

		if ( s == "ndXcacdd" )  ndXcacdd    = ToDouble ( xml->GetText() );

		if ( s == "ndXcacbd" )  ndXcacbd    = ToDouble ( xml->GetText() );

		if ( s == "ndYvp" )  ndYvp       = ToDouble ( xml->GetText() );

		if ( s == "ndYur" )  ndYur       = ToDouble ( xml->GetText() );

		if ( s == "ndYuv" )  ndYuv       = ToDouble ( xml->GetText() );

		if ( s == "ndYvav" )  ndYvav      = ToDouble ( xml->GetText() );

		if ( s == "ndYcacdad" )  ndYcacdad   = ToDouble ( xml->GetText() );

		if ( s == "ndYcacbabad" )  ndYcacbabad = ToDouble ( xml->GetText() );

		if ( s == "ndYT" )  ndYT        = ToDouble ( xml->GetText() );

		if ( s == "ndNrp" )  ndNrp       = ToDouble ( xml->GetText() );

		if ( s == "ndNur" )  ndNur       = ToDouble ( xml->GetText() );

		if ( s == "ndNuv" )  ndNuv       = ToDouble ( xml->GetText() );

		if ( s == "ndNanr" )  ndNanr      = ToDouble ( xml->GetText() );

		if ( s == "ndNcacd" )  ndNcacd     = ToDouble ( xml->GetText() );

		if ( s == "ndNcacbabad" )  ndNcacbabad = ToDouble ( xml->GetText() );

		if ( s == "ndNT" )  ndNT        = ToDouble ( xml->GetText() );

		if ( s == "ndTuu" )  ndTuu       = ToDouble ( xml->GetText() );

		if ( s == "ndTun" )  ndTun       = ToDouble ( xml->GetText() );

		if ( s == "ndTnan" )  ndTnan      = ToDouble ( xml->GetText() );

		if ( s == "ndCun" )  ndCun       = ToDouble ( xml->GetText() );

		if ( s == "ndCnn" )  ndCnn       = ToDouble ( xml->GetText() );

		if ( s == "tt" )  tt          = ToDouble ( xml->GetText() );

		xml = n->GetNextChild ( "DOUBLE" );
	}
}

void Model2v1::init() {
	g    = 9.80665; 			//Acceleration due to gravity (m/s^2)
	rho  = 1000.; 			//denisty of water (kg/m^3)


	d     = ( dfor + daft ) / 2; 		//Draft (m)
	m     = rho * Cb * LWL * d * B; 	//Mass (kg)
	Izz   = m * pow ( Rg,2 ); 		//Moment of inertia


	Xup       = ndXup       * m;
	Xuu       = ndXuu       * m/L;
	Xvr       = ndXvr       * m;
	Xvv       = ndXvv       * m/L;
	Xcacdd    = ndXcacdd    * m/L;
	Xcacbd    = ndXcacbd    * m/L;
	Yvp       = ndYvp       * m;
	Yur       = ndYur       * m;
	Yuv       = ndYuv       * m/L;
	Yvav      = ndYvav      * m/L;
	Ycacdad   = ndYcacdad   * m/L;
	Ycacbabad = ndYcacbabad * m/L;
	YT        = ndYT        * 1;
	Nrp       = ndNrp       * m*L*L;
	Nur       = ndNur       * m*L;
	Nuv       = ndNuv       * m;
	Nanr      = ndNanr      * m*L;
	Ncacd     = ndNcacd     * m;
	Ncacbabad = ndNcacbabad * m;
	NT        = ndNT        * L;
	Tuu       = ndTuu       * m/L;
	Tun       = ndTun       * m/L;
	Tnan      = ndTnan      * m/L;
	Cun       = ndCun       * sqrt ( L/g ) / ( 2.*M_PI );
	Cnn       = ndCnn       * sqrt ( L/g ) *L/ ( 2.*M_PI*2*M_PI );
}



