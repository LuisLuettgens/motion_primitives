#include "Model1v2.h"
using namespace std;



Model1v2::Model1v2 ( string xmlparams ) :BaseModel ( xmlparams ) {

	loadParams ( xmlparams );
	init();

	delta. resize ( 2 );
	n.     resize ( 2 );
	deltaC.resize ( 2 );
	EOT.   resize ( 2 );
	PP.    resize ( 2 );
	kT.    resize ( 1 );
}

void Model1v2::ode ( double *dx, const double *x, const double *ctrl, const double *p, const double *q ) {
	x0        = x[0];
	y0        = x[1];
	psi       = x[2];
	u         = x[3];
	v         = x[4];
	r         = x[5];
	delta[0]  = x[6];
	delta[1]  = x[7];
	n[0]      = x[8];
	n[1]      = x[9];

	deltaC[0] = 0;
	deltaC[1] = 0;
	EOT[0]    = 0;
	EOT[1]    = 0;
	PP[0]  = 0;
	PP[1]  = 0;
	kT[0]  = 0;

	DFT  = 1 * 463./900.;//q[0];
	SET  = M_PI/2;//q[1];
	TWS  = 15 * 463./900.;//q[2];
	TWA  = q[3];
	TWD  = q[4];
	AWS  = q[5];
	AWA  = q[6];
	AWD  = q[7];
	h    = q[8];

	dx[0] = u * cos ( psi ) - v * sin ( psi ) + SET*cos ( DFT );
	dx[1] = u * sin ( psi ) + v * cos ( psi ) + SET*sin ( DFT );
	dx[2] = r;

	//Ideal fluid effects
	XI = Xvr*v*r + Xrr*r*r + Xvv*v*v;
	YI = Xup*u*r;
	NI = Nvp*u*r + ( Yvp - Xup ) *u*v;


	//Hull lifting effects
	vs = k*r*l*sign ( u );
	sqroot = sqrt ( u*u + ( vs-v ) * ( vs-v ) );

	if ( sqroot < 1e-15 ) {
		XHL = 0;
		YHL = 0;
		NHL = 0;
	} else {
		XHL =  rho/2*L*  d * u* ( b1p*vs-b1*v ) /sqroot * ( ( vs-v ) - u*u* ( b2p*vs-b2*v ) / ( sqroot*sqroot ) );
		YHL =  rho/2*L*  d * u*u* ( b1p*vs-b1*v ) /sqroot * ( 1 + ( b2p*vs-b2*v ) * ( vs-v ) / ( sqroot*sqroot ) );
		NHL = -rho/2*L*L*d * u*abs ( u ) * ( b3p*vs-b3*v ) /sqroot * ( 1 + ( b2p*vs-b2*v ) * ( vs-v ) / ( sqroot*sqroot ) );
	}


	//Hull longitudinal resistance
	XHR = - ( XHRu * u + XHRuau * u * abs ( u ) + XHRuuu * u*u*u );

	//Hull cross flow
	if ( v*v >= r*r * l*l ) {
		YHC = -rho/2.*L*  d * ( ( a0+a8/9. ) *pow ( v,2 ) +
								( a0/3.+a8/11. ) *pow ( r*l,2 ) +
								2.* ( a7/9. + a9/11. ) *v*r*l ) * sign ( v );

		NHC = -rho/2.*L*L*d * ( ( a7/18.+a9/22. ) *pow ( v,2 ) +
								( a7/22.+a9/26. ) *pow ( r*l,2 ) +
								( a0/3. + a8/11. ) *v*r*l ) * sign ( v );

	} else {
		YHC = -rho/2.*L*  d * ( a0* ( r*l + v*v/ ( 3.*r*l ) ) *v +
								a7/2.* ( v*v/4. + r*r*l*l/5. - v*v/180.*pow ( v/ ( r*l ),8 ) ) +
								a8/5.* ( r*l + v/99.*pow ( v/ ( r*l ),9 ) ) *v +
								a9/2.* ( v*v/5. + r*r*l*l/6. - v*v/330.*pow ( v/ ( r*l ),10 ) ) ) * sign ( r );

		NHC = -rho/2.*L*L*d * ( a0/4.* ( v*v + r*r*l*l/2. - pow ( v,4 ) / ( 6.*r*r*l*l ) ) +
								a7/10.* ( r*l + v/99.*pow ( v/ ( r*l ),9 ) ) *v +
								a8/4.* ( v*v/5. + r*r*l*l/6. - v*v/330.*pow ( v/ ( r*l ),10 ) ) +
								a9/12.* ( r*l + v/143.*pow ( v/ ( r*l ),11 ) ) *v ) * sign ( r );
	}


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

	//Propellor1
	cP1 = 0.7 * M_PI * n[0] * DP;
	uP1 = ( 1 - omP ) * u;

	if ( abs ( uP1 ) < 1e-15 and abs ( cP1 ) < 1e-15 ) {
		T1 = 0;
		Q1 = 0;
	} else {
		eps1 = atan2 ( uP1,cP1 );

		if ( n[0] < 0 ) {
			YPT1 =  YPTm;
			NPT1 =  NPTm;
		} else {
			YPT1 = YPTp;
			NPT1 = NPTp;
		};

		CT1 = 0;

		for ( int j=0; j<21; j++ ) CT1 += AKT[j]*cos ( j*eps1 ) + BKT[j]*sin ( j*eps1 );

		T1 = rho/2.*AP* ( uP1*uP1 + cP1*cP1 ) *CT1;
		CQ1 = 0;

		for ( int j=0; j<21; j++ ) CQ1 += AKQ[j]*cos ( j*eps1 ) + BKQ[j]*sin ( j*eps1 );

		Q1 = rho/2.*AP*DP* ( uP1*uP1 + cP1*cP1 ) *CQ1;
	}

	XP1 = ( 1 - tP ) * T1;
	YP1 = RotS1 * YPT1 * T1;
	NP1 = RotS1 * NPT1 * T1 + yP*T1;

	//Propellor2
	cP2 = 0.7 * M_PI * n[1] * DP;
	uP2 = ( 1 - omP ) * u;

	if ( abs ( uP2 ) < 1e-15 and abs ( cP2 ) < 1e-15 ) {
		T2 = 0;
		Q2 = 0;
	} else {
		eps2 = atan2 ( uP2,cP2 );

		if ( n[1] < 0 ) {
			YPT2 = YPTm;
			NPT2 = NPTm;
		} else {
			YPT2 = YPTp;
			NPT2 = NPTp;
		};

		CT2 = 0;

		for ( int j=0; j<21; j++ ) CT2 += AKT[j]*cos ( j*eps2 ) + BKT[j]*sin ( j*eps2 );

		T2 = rho/2.*AP* ( uP2*uP2 + cP2*cP2 ) *CT2;
		CQ2 = 0;

		for ( int j=0; j<21; j++ ) CQ2 += AKQ[j]*cos ( j*eps2 ) + BKQ[j]*sin ( j*eps2 );

		Q2 = rho/2.*AP*DP* ( uP2*uP2 + cP2*cP2 ) *CQ2;
	}

	XP2 = ( 1 - tP ) * T2;
	YP2 = RotS2 * YPT2 * T2;
	NP2 = RotS2 * NPT2 * T2 - yP*T2;



	//Rudder
	lift0 = 0;
	lift5 = 0;
	drag0 = 0;
	drag5 = 0;
	de0 = 0;
	de5 = M_PI;
	lia1 = ( lift1-lift0 ) / ( de1-de0 ) - lia2* ( de0+de1 ) - lia3* ( de0*de0 + de0*de1 + de1*de1 );
	lia0 = lift0 - lia1*de0 - lia2*de0*de0 - lia3*de0*de0*de0;
	lib1 = ( lift5-lift4 ) / ( de5-de4 ) - lib2* ( de4+de5 ) - lib3* ( de4*de4 + de4*de5 + de5*de5 );
	lib0 = lift4 - lib1*de4 - lib2*de4*de4 - lib3*de4*de4*de4;
	dra1 = ( drag1-drag0 ) / ( de1-de0 ) - dra2* ( de0+de1 ) - dra3* ( de0*de0 + de0*de1 + de1*de1 );
	dra0 = drag0 - dra1*de0 - dra2*de0*de0 - dra3*de0*de0*de0;
	drb1 = ( drag3-drag2 ) / ( de3-de2 ) - dra2* ( de2+de3 ) - dra3* ( de2*de2 + de2*de3 + de3*de3 );
	drb0 = drag2 - drb1*de2 - drb2*de2*de2 - drb3*de2*de2*de2;
	drc1 = ( drag5-drag4 ) / ( de5-de4 ) - dra2* ( de4+de5 ) - dra3* ( de4*de4 + de4*de5 + de5*de5 );
	drc0 = drag4 - drc1*de4 - drc2*de4*de4 - drc3*de4*de4*de4;

	//Rudder1
	uR1 = ( 1 - omR ) * u;
	vR1 = ( v + r * xR ) * kHR;

	uRP1 = uP1 + sign ( T1 ) * ( ( kPR - 1./2 ) * sign ( T1 ) + 1./2 ) * ( sqrt ( pow ( uP1,2 ) + 2*abs ( T1 ) / ( rho*AP ) ) - sign ( u ) *uP1 );

	if ( uRP1*uP1 >= 0 ) {
		uRmean1 = sign ( uP1 ) *sqrt ( eta*uRP1*uRP1 + ( 1-eta ) *uP1*uP1 );
	} else {
		uRmean1 = sign ( uRP1* ( eta*uRP1*uRP1 - ( 1-eta ) *uP1*uP1 ) ) *sqrt ( abs ( eta*uRP1*uRP1 - ( 1-eta ) *uP1*uP1 ) );
	}

	betaR1 = atan2 ( -vR1, uRmean1 );
	deltae1 = delta[0] + betaR1;
	deltae1 = atan2 ( sin ( deltae1 ),cos ( deltae1 ) );

	if ( abs ( deltae1 ) >= de0 && abs ( deltae1 ) < de1 ) {
		CLR01 = lia0 + lia1*abs ( deltae1 ) + lia2*abs ( deltae1 ) *abs ( deltae1 ) + lia3*abs ( deltae1 ) *abs ( deltae1 ) *abs ( deltae1 );
		CDR01 = dra0 + dra1*abs ( deltae1 ) + dra2*abs ( deltae1 ) *abs ( deltae1 ) + dra3*abs ( deltae1 ) *abs ( deltae1 ) *abs ( deltae1 );
	} else if ( abs ( deltae1 ) >= de1 && abs ( deltae1 ) < de2 ) {
		CLR01 = ( lift2-lift1 ) / ( de2-de1 ) *abs ( deltae1 ) + ( de2*lift1-de1*lift2 ) / ( de2-de1 );
		CDR01 = ( drag2-drag1 ) / ( de2-de1 ) *abs ( deltae1 ) + ( de2*drag1-de1*drag2 ) / ( de2-de1 );
	} else if ( abs ( deltae1 ) >= de2 && abs ( deltae1 ) < de3 ) {
		CLR01 = ( lift3-lift2 ) / ( de3-de2 ) *abs ( deltae1 ) + ( de3*lift2-de2*lift3 ) / ( de3-de2 );
		CDR01 = drb0 + drb1*abs ( deltae1 ) + drb2*abs ( deltae1 ) *abs ( deltae1 ) + drb3*abs ( deltae1 ) *abs ( deltae1 ) *abs ( deltae1 );
	} else if ( abs ( deltae1 ) >= de3 && abs ( deltae1 ) < de4 ) {
		CLR01 = ( lift4-lift3 ) / ( de4-de3 ) *abs ( deltae1 ) + ( de4*lift3-de3*lift4 ) / ( de4-de3 );
		CDR01 = ( drag4-drag3 ) / ( de4-de3 ) *abs ( deltae1 ) + ( de4*drag3-de3*drag4 ) / ( de4-de3 );
	} else if ( abs ( deltae1 ) >= de4 && abs ( deltae1 ) <= de5 ) {
		CLR01 = lib0 + lib1*abs ( deltae1 ) + lib2*abs ( deltae1 ) *abs ( deltae1 ) + lib3*abs ( deltae1 ) *abs ( deltae1 ) *abs ( deltae1 );
		CDR01 = drc0 + drc1*abs ( deltae1 ) + drc2*abs ( deltae1 ) *abs ( deltae1 ) + drc3*abs ( deltae1 ) *abs ( deltae1 ) *abs ( deltae1 );
	}

	if ( deltae1 < 0 ) CLR01 = -CLR01;

	CLR01 *= kCLR;
	CDR01 *= kCDR;
	CLR1 = ( 1 + kLR/ ( 1+ ( uR1-uRmean1 ) * ( uR1-uRmean1 ) ) ) *CLR01;
	CDR1 = ( 1 + kDR/ ( 1+ ( uR1-uRmean1 ) * ( uR1-uRmean1 ) ) ) *CDR01;

	XR1 = rho / 2 * AR * ( pow ( uRmean1,2 ) + pow ( vR1,2 ) ) * ( CLR1 * sin ( betaR1 ) - CDR1 * cos ( betaR1 ) );
	YR1 = rho / 2 * AR * ( pow ( uRmean1,2 ) + pow ( vR1,2 ) ) * ( CLR1 * cos ( betaR1 ) + CDR1 * sin ( betaR1 ) );
	NR1 = ( 1 - kNR/ ( 1+ ( uR1-uRmean1 ) * ( uR1-uRmean1 ) ) ) *YR1 * xR;


	//Rudder2
	uR2 = ( 1 - omR ) * u;
	vR2 = ( v + r * xR ) * kHR;

	uRP2 = uP2 + sign ( T2 ) * ( ( kPR - 1./2 ) * sign ( T2 ) + 1./2 ) * ( sqrt ( pow ( uP2,2 ) + 2*abs ( T2 ) / ( rho*AP ) ) - sign ( u ) *uP2 );

	if ( uRP2*uP2 >= 0 ) {
		uRmean2 = sign ( uP2 ) *sqrt ( eta*uRP2*uRP2 + ( 1-eta ) *uP2*uP2 );
	} else {
		uRmean2 = sign ( uRP2* ( eta*uRP2*uRP2 - ( 1-eta ) *uP2*uP2 ) ) *sqrt ( abs ( eta*uRP2*uRP2 - ( 1-eta ) *uP2*uP2 ) );
	}

	betaR2 = atan2 ( -vR2, uRmean2 );
	deltae2 = delta[1] + betaR2;
	deltae2 = atan2 ( sin ( deltae2 ),cos ( deltae2 ) );

	if ( abs ( deltae2 ) >= de0 && abs ( deltae2 ) < de1 ) {
		CLR02 = lia0 + lia1*abs ( deltae2 ) + lia2*abs ( deltae2 ) *abs ( deltae2 ) + lia3*abs ( deltae2 ) *abs ( deltae2 ) *abs ( deltae2 );
		CDR02 = dra0 + dra1*abs ( deltae2 ) + dra2*abs ( deltae2 ) *abs ( deltae2 ) + dra3*abs ( deltae2 ) *abs ( deltae2 ) *abs ( deltae2 );
	} else if ( abs ( deltae2 ) >= de1 && abs ( deltae2 ) < de2 ) {
		CLR02 = ( lift2-lift1 ) / ( de2-de1 ) *abs ( deltae2 ) + ( de2*lift1-de1*lift2 ) / ( de2-de1 );
		CDR02 = ( drag2-drag1 ) / ( de2-de1 ) *abs ( deltae2 ) + ( de2*drag1-de1*drag2 ) / ( de2-de1 );
	} else if ( abs ( deltae2 ) >= de2 && abs ( deltae2 ) < de3 ) {
		CLR02 = ( lift3-lift2 ) / ( de3-de2 ) *abs ( deltae2 ) + ( de3*lift2-de2*lift3 ) / ( de3-de2 );
		CDR02 = drb0 + drb1*abs ( deltae2 ) + drb2*abs ( deltae2 ) *abs ( deltae2 ) + drb3*abs ( deltae2 ) *abs ( deltae2 ) *abs ( deltae2 );
	} else if ( abs ( deltae2 ) >= de3 && abs ( deltae2 ) < de4 ) {
		CLR02 = ( lift4-lift3 ) / ( de4-de3 ) *abs ( deltae2 ) + ( de4*lift3-de3*lift4 ) / ( de4-de3 );
		CDR02 = ( drag4-drag3 ) / ( de4-de3 ) *abs ( deltae2 ) + ( de4*drag3-de3*drag4 ) / ( de4-de3 );
	} else if ( abs ( deltae2 ) >= de4 && abs ( deltae2 ) <= de5 ) {
		CLR02 = lib0 + lib1*abs ( deltae2 ) + lib2*abs ( deltae2 ) *abs ( deltae2 ) + lib3*abs ( deltae2 ) *abs ( deltae2 ) *abs ( deltae2 );
		CDR02 = drc0 + drc1*abs ( deltae2 ) + drc2*abs ( deltae2 ) *abs ( deltae2 ) + drc3*abs ( deltae2 ) *abs ( deltae2 ) *abs ( deltae2 );
	}

	if ( deltae2 < 0 ) CLR02 = -CLR02;

	CLR02 *= kCLR;
	CDR02 *= kCDR;
	CLR2 = ( 1 + kLR/ ( 1+ ( uR2-uRmean2 ) * ( uR2-uRmean2 ) ) ) *CLR02;
	CDR2 = ( 1 + kDR/ ( 1+ ( uR2-uRmean2 ) * ( uR2-uRmean2 ) ) ) *CDR02;

	XR2 = rho / 2 * AR * ( pow ( uRmean2,2 ) + pow ( vR2,2 ) ) * ( CLR2 * sin ( betaR2 ) - CDR2 * cos ( betaR2 ) );
	YR2 = rho / 2 * AR * ( pow ( uRmean2,2 ) + pow ( vR2,2 ) ) * ( CLR2 * cos ( betaR2 ) + CDR2 * sin ( betaR2 ) );
	NR2 = ( 1 - kNR/ ( 1+ ( uR2-uRmean2 ) * ( uR2-uRmean2 ) ) ) *YR2 * xR;

	//Thruster
	kh = min ( 1.0, ah + bh*hT/RT );

	wjet = sqrt ( abs ( TT ) / ( rho*AT ) );

	uTrel = u/wjet;

	aTY1 = ( kYm-1 ) / ( uTYrel1-uTYrel2 )-aTY2* ( uTYrel1+uTYrel2 ) - aTY3* ( uTYrel1*uTYrel1 + uTYrel1*uTYrel2 + uTYrel2*uTYrel2 );
	aTY0 = 1 - aTY1*uTYrel2 - aTY2*uTYrel2*uTYrel2 - aTY3*uTYrel2*uTYrel2*uTYrel2;
	bTY1 = ( kYp-1 ) / ( uTYrel4-uTYrel3 )-bTY2* ( uTYrel3+uTYrel4 ) - bTY3* ( uTYrel3*uTYrel3 + uTYrel3*uTYrel4 + uTYrel4*uTYrel4 );
	bTY0 = 1 - bTY1*uTYrel3 - bTY2*uTYrel3*uTYrel3 - bTY3*uTYrel3*uTYrel3*uTYrel3;

	aTN1 = ( kNm-1 ) / ( uTNrel1-uTNrel2 )-aTN2* ( uTNrel1+uTNrel2 ) - aTN3* ( uTNrel1*uTNrel1 + uTNrel1*uTNrel2 + uTNrel2*uTNrel2 );
	aTN0 = 1 - aTN1*uTNrel2 - aTN2*uTNrel2*uTNrel2 - aTN3*uTNrel2*uTNrel2*uTNrel2;
	bTN1 = ( kNp-1 ) / ( uTNrel4-uTNrel3 )-bTN2* ( uTNrel3+uTNrel4 ) - bTN3* ( uTNrel3*uTNrel3 + uTNrel3*uTNrel4 + uTNrel4*uTNrel4 );
	bTN0 = 1 - bTN1*uTNrel3 - bTN2*uTNrel3*uTNrel3 - bTN3*uTNrel3*uTNrel3*uTNrel3;

	if ( uTrel < uTYrel1 ) {
		kY = kYm;
	} else if ( uTYrel1 <= uTrel && uTrel < uTYrel2 ) {
		kY = aTY0 + aTY1*uTrel + aTY2*uTrel*uTrel + aTY3*uTrel*uTrel*uTrel;
	} else if ( uTYrel2 <= uTrel && uTrel < uTYrel3 ) {
		kY = 1;
	} else if ( uTYrel3 <= uTrel && uTrel < uTYrel4 ) {
		kY = bTY0 + bTY1*uTrel + bTY2*uTrel*uTrel + bTY3*uTrel*uTrel*uTrel;
	} else if ( uTYrel4 <= uTrel ) {
		kY = kYp;
	}

	if ( uTrel < uTNrel1 ) {
		kN = kNm;
	} else if ( uTNrel1 <= uTrel && uTrel < uTNrel2 ) {
		kN = aTN0 + aTN1*uTrel + aTN2*uTrel*uTrel + aTN3*uTrel*uTrel*uTrel;
	} else if ( uTNrel2 <= uTrel && uTrel < uTNrel3 ) {
		kN = 1;
	} else if ( uTNrel3 <= uTrel && uTrel < uTNrel4 ) {
		kN = bTN0 + bTN1*uTrel + bTN2*uTrel*uTrel + bTN3*uTrel*uTrel*uTrel;
	} else if ( uTNrel4 <= uTrel ) {
		kN = kNp;
	}

	YT = kh*kY*kT[0]*TTmax;
	NT = kh*kN*kT[0]*TTmax*xT;


	//Engine
	IEP = 1.e+6;
	QEmax1 = 4.1e+7;
	QEmax2 = 4.1e+7;

	double a = 0.075;
	double b = 0.250;
	double Af = 2.5;
	double Bf = 1.0;
	double Ab = 1.0;
	double Bb = 0.6;

	if ( EOT[0] >= 0 ) QE1 = Af* ( EOT[0]/1.-a ) / ( 1.-a ) * ( 1.-n[0]/1.5 ) + Bf* ( EOT[0]/1.-b ) / ( 1.-b ) *n[0]/1.5 + 15./74;

	if ( EOT[0] <  0 ) QE1 = Ab* ( EOT[0]/1.+a ) / ( 1.-a ) * ( 1.+n[0]/1.5 ) - Bb* ( EOT[0]/1.+b ) / ( 1.-b ) *n[0]/1.5 - 15./74;

	if ( EOT[1] >= 0 ) QE2 = Af* ( EOT[1]/1.-a ) / ( 1.-a ) * ( 1.-n[1]/1.5 ) + Bf* ( EOT[1]/1.-b ) / ( 1.-b ) *n[1]/1.5 + 15./74;

	if ( EOT[1] <  0 ) QE2 = Ab* ( EOT[1]/1.+a ) / ( 1.-a ) * ( 1.+n[1]/1.5 ) - Bb* ( EOT[1]/1.+b ) / ( 1.-b ) *n[1]/1.5 - 15./74;

	QE1 *= QEmax1;
	QE2 *= QEmax2;




	dx[3] = XI + XHL + XHR + XP1 + XP2 + XR1 + XR2      + XA  +  m*v*r + m*xG*r*r;
	dv    = YI + YHL + YHC + YP1 + YP2 + YR1 + YR2 + YT + YA  -  m*u*r           ;
	dr    = NI + NHL + NHC + NP1 + NP2 + NR1 + NR2 + NT + NA  -  m*xG*u*r        ;

	dx[3] = 1./m11 * dx[3];
	dx[4] = 1./det * ( m33*dv - m23*dr );
	dx[5] = 1./det * ( -m32*dv + m22*dr );
	dx[6] = 1./LagR* ( deltaC[0] - delta[0] );
	dx[7] = 1./LagR* ( deltaC[1] - delta[1] );
	dx[8] = 1./ ( 2*M_PI*IEP ) * ( QE1 - Q1 ) *1.e-1;
	dx[9] = 1./ ( 2*M_PI*IEP ) * ( QE2 - Q2 ) *1.e-1;

//       cout << uP1 << " " << n[0] << " " << eps1*180/M_PI << " " << QE1 - Q1 << " " << dx[8] << endl;

}

void Model1v2::ode_structure ( tw::DiffStructure &s, int N ) {
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
	s ( N+3, N+6 );
	s ( N+3, N+7 );
	s ( N+3, N+8 );
	s ( N+3, N+9 );
	s ( N+4, N+3 );
	s ( N+4, N+4 );
	s ( N+4, N+5 );
	s ( N+4, N+6 );
	s ( N+4, N+7 );
	s ( N+4, N+8 );
	s ( N+4, N+9 );
	s ( N+5, N+3 );
	s ( N+5, N+4 );
	s ( N+5, N+5 );
	s ( N+5, N+6 );
	s ( N+5, N+7 );
	s ( N+5, N+8 );
	s ( N+5, N+9 );

	int nbrS = get_nbrStates();
	int nbrC = get_nbrControls();

	for ( int i=0; i<nbrC; i++ ) {
		s ( N+3, N+nbrS+i );
		s ( N+4, N+nbrS+i );
		s ( N+5, N+nbrS+i );
	}

	s ( N+6, N+6 );
	s ( N+6, N+nbrS+0 );
	s ( N+7, N+7 );
	s ( N+7, N+nbrS+1 );

	s ( N+8, N+3 );
	s ( N+8, N+8 );
	s ( N+8, N+nbrS+2 );
	s ( N+9, N+3 );
	s ( N+9, N+9 );
	s ( N+9, N+nbrS+3 );
}

void Model1v2::loadParams ( string xmlparams ) {
	XMLParser p;
	XMLNode *n = p.Parse ( xmlparams )->GetFirstChild ( "PARAMS" );
	XMLNode *xml = n->GetFirstChild ( "DOUBLE" );

	while ( xml ) {
		string s ( xml->GetAttribute ( "name" ) );

		if ( s == "g" ) g         = ToDouble ( xml->GetText() );

		if ( s == "rho" ) rho       = ToDouble ( xml->GetText() );

		if ( s == "rhoA" ) rhoA      = ToDouble ( xml->GetText() );

		//Hull
		if ( s == "L" )  L        = ToDouble ( xml->GetText() );

		if ( s == "LWL" )  LWL      = ToDouble ( xml->GetText() );

		if ( s == "LOA" )  LOA      = ToDouble ( xml->GetText() );

		if ( s == "AF" )  AF       = ToDouble ( xml->GetText() );

		if ( s == "AL" )  AL       = ToDouble ( xml->GetText() );

		if ( s == "BOA" )  BOA      = ToDouble ( xml->GetText() );

		if ( s == "dfor" )  dfor     = ToDouble ( xml->GetText() );

		if ( s == "daft" )  daft     = ToDouble ( xml->GetText() );

		if ( s == "Cb" )  Cb       = ToDouble ( xml->GetText() );

		if ( s == "xG" )  xG       = ToDouble ( xml->GetText() );

		if ( s == "Rg" )  Rg       = ToDouble ( xml->GetText() );

		if ( s == "ndXup" )  ndXup    = ToDouble ( xml->GetText() );

		if ( s == "ndXrr" )  ndXrr    = ToDouble ( xml->GetText() );

		if ( s == "ndYvp" )  ndYvp    = ToDouble ( xml->GetText() );

		if ( s == "ndNrp" )  ndNrp    = ToDouble ( xml->GetText() );

		if ( s == "ndXvr" )  ndXvr    = ToDouble ( xml->GetText() );

		if ( s == "ndXvv" )  ndXvv    = ToDouble ( xml->GetText() );

		if ( s == "ndYrp" )  ndYrp    = ToDouble ( xml->GetText() );

		if ( s == "ndNvp" )  ndNvp    = ToDouble ( xml->GetText() );

		if ( s == "b1" )  b1       = ToDouble ( xml->GetText() );

		if ( s == "b2" )  b2       = ToDouble ( xml->GetText() );

		if ( s == "b3" )  b3       = ToDouble ( xml->GetText() );

		if ( s == "b1p" )  b1p      = ToDouble ( xml->GetText() );

		if ( s == "b2p" )  b2p      = ToDouble ( xml->GetText() );

		if ( s == "b3p" )  b3p      = ToDouble ( xml->GetText() );

		if ( s == "k" )  k        = ToDouble ( xml->GetText() );

		if ( s == "a0" )  a0       = ToDouble ( xml->GetText() );

		if ( s == "a7" )  a7       = ToDouble ( xml->GetText() );

		if ( s == "a8" )  a8       = ToDouble ( xml->GetText() );

		if ( s == "a9" )  a9       = ToDouble ( xml->GetText() );

		if ( s == "ndXHRu" ) ndXHRu   = ToDouble ( xml->GetText() );

		if ( s == "ndXHRuau" ) ndXHRuau = ToDouble ( xml->GetText() );

		if ( s == "ndXHRuuu" ) ndXHRuuu = ToDouble ( xml->GetText() );

		//Propeller
		if ( s == "DP" )  DP       = ToDouble ( xml->GetText() );

		if ( s == "tP" )  tP       = ToDouble ( xml->GetText() );

		if ( s == "omP" )  omP      = ToDouble ( xml->GetText() );

		if ( s == "prT0" )  prT[0]   = ToDouble ( xml->GetText() );

		if ( s == "prT1" )  prT[1]   = ToDouble ( xml->GetText() );

		if ( s == "prT2" )  prT[2]   = ToDouble ( xml->GetText() );

		if ( s == "prT3" )  prT[3]   = ToDouble ( xml->GetText() );

		if ( s == "prT4" )  prT[4]   = ToDouble ( xml->GetText() );

		if ( s == "prT5" )  prT[5]   = ToDouble ( xml->GetText() );

		if ( s == "prT6" )  prT[6]   = ToDouble ( xml->GetText() );

		if ( s == "prT7" )  prT[7]   = ToDouble ( xml->GetText() );

		if ( s == "prT8" )  prT[8]   = ToDouble ( xml->GetText() );

		if ( s == "prT9" )  prT[9]   = ToDouble ( xml->GetText() );

		if ( s == "prT10" )  prT[10]  = ToDouble ( xml->GetText() );

		if ( s == "ndxP" )  ndxP     = ToDouble ( xml->GetText() );

		if ( s == "ndyP" )  ndyP     = ToDouble ( xml->GetText() );

		if ( s == "ndYPTp" )  ndYPTp   = ToDouble ( xml->GetText() );

		if ( s == "ndYPTm" )  ndYPTm   = ToDouble ( xml->GetText() );

		if ( s == "ndNPTp" )  ndNPTp   = ToDouble ( xml->GetText() );

		if ( s == "ndNPTm" )  ndNPTm   = ToDouble ( xml->GetText() );

		if ( s == "RotS1" )  RotS1    = ToDouble ( xml->GetText() );

		if ( s == "RotS2" )  RotS2    = ToDouble ( xml->GetText() );

		//Rudder
		if ( s == "LagR" )  LagR     = ToDouble ( xml->GetText() );

		if ( s == "AR" )  AR       = ToDouble ( xml->GetText() );

		if ( s == "ChordR" )  ChordR   = ToDouble ( xml->GetText() );

		if ( s == "omR" )  omR      = ToDouble ( xml->GetText() );

		if ( s == "kHR" )  kHR      = ToDouble ( xml->GetText() );

		if ( s == "kPR" )  kPR      = ToDouble ( xml->GetText() );

		if ( s == "ndxR" )  ndxR     = ToDouble ( xml->GetText() );

		if ( s == "kLR" )  kLR      = ToDouble ( xml->GetText() );

		if ( s == "kDR" )  kDR      = ToDouble ( xml->GetText() );

		if ( s == "kNR" )  kNR      = ToDouble ( xml->GetText() );

		if ( s == "kCLR" )  kCLR     = ToDouble ( xml->GetText() );

		if ( s == "kCDR" )  kCDR     = ToDouble ( xml->GetText() );

		if ( s == "lift1" ) lift1     = ToDouble ( xml->GetText() );

		if ( s == "lift2" ) lift2     = ToDouble ( xml->GetText() );

		if ( s == "lift3" ) lift3     = ToDouble ( xml->GetText() );

		if ( s == "lift4" ) lift4     = ToDouble ( xml->GetText() );

		if ( s == "drag1" ) drag1     = ToDouble ( xml->GetText() );

		if ( s == "drag2" ) drag2     = ToDouble ( xml->GetText() );

		if ( s == "drag3" ) drag3     = ToDouble ( xml->GetText() );

		if ( s == "drag4" ) drag4     = ToDouble ( xml->GetText() );

		if ( s == "lia2" ) lia2      = ToDouble ( xml->GetText() );

		if ( s == "lia3" ) lia3      = ToDouble ( xml->GetText() );

		if ( s == "lib2" ) lib2      = ToDouble ( xml->GetText() );

		if ( s == "lib3" ) lib3      = ToDouble ( xml->GetText() );

		if ( s == "dra2" ) dra2      = ToDouble ( xml->GetText() );

		if ( s == "dra3" ) dra3      = ToDouble ( xml->GetText() );

		if ( s == "drb2" ) drb2      = ToDouble ( xml->GetText() );

		if ( s == "drb3" ) drb3      = ToDouble ( xml->GetText() );

		if ( s == "drc2" ) drc2      = ToDouble ( xml->GetText() );

		if ( s == "drc3" ) drc3      = ToDouble ( xml->GetText() );

		if ( s == "de1" ) de1       = ToDouble ( xml->GetText() );

		if ( s == "de2" ) de2       = ToDouble ( xml->GetText() );

		if ( s == "de3" ) de3       = ToDouble ( xml->GetText() );

		if ( s == "de4" ) de4       = ToDouble ( xml->GetText() );

		//Wind
		if ( s == "CX0" )  CX0      = ToDouble ( xml->GetText() );

		if ( s == "CY0" )  CY0      = ToDouble ( xml->GetText() );

		if ( s == "xA0" )  xA0      = ToDouble ( xml->GetText() );

		//Thruster
		if ( s == "ah" )  ah       = ToDouble ( xml->GetText() );

		if ( s == "bh" )  bh       = ToDouble ( xml->GetText() );

		if ( s == "hT" )  hT       = ToDouble ( xml->GetText() );

		if ( s == "RT" )  RT       = ToDouble ( xml->GetText() );

		if ( s == "AT" )  AT       = ToDouble ( xml->GetText() );

		if ( s == "xT" )  xT       = ToDouble ( xml->GetText() );

		if ( s == "aTY2" )  aTY2     = ToDouble ( xml->GetText() );

		if ( s == "aTY3" )  aTY3     = ToDouble ( xml->GetText() );

		if ( s == "bTY2" )  bTY2     = ToDouble ( xml->GetText() );

		if ( s == "bTY3" )  bTY3     = ToDouble ( xml->GetText() );

		if ( s == "aTN2" )  aTN2     = ToDouble ( xml->GetText() );

		if ( s == "aTN3" )  aTN3     = ToDouble ( xml->GetText() );

		if ( s == "bTN2" )  bTN2     = ToDouble ( xml->GetText() );

		if ( s == "bTN3" )  bTN3     = ToDouble ( xml->GetText() );

		if ( s == "kYm" )  kYm      = ToDouble ( xml->GetText() );

		if ( s == "kYp" )  kYp      = ToDouble ( xml->GetText() );

		if ( s == "kNm" )  kNm      = ToDouble ( xml->GetText() );

		if ( s == "kNp" )  kNp      = ToDouble ( xml->GetText() );

		if ( s == "uTYrel1" )  uTYrel1  = ToDouble ( xml->GetText() );

		if ( s == "uTYrel2" )  uTYrel2  = ToDouble ( xml->GetText() );

		if ( s == "uTYrel3" )  uTYrel3  = ToDouble ( xml->GetText() );

		if ( s == "uTYrel4" )  uTYrel4  = ToDouble ( xml->GetText() );

		if ( s == "uTNrel1" )  uTNrel1  = ToDouble ( xml->GetText() );

		if ( s == "uTNrel2" )  uTNrel2  = ToDouble ( xml->GetText() );

		if ( s == "uTNrel3" )  uTNrel3  = ToDouble ( xml->GetText() );

		if ( s == "uTNrel4" )  uTNrel4  = ToDouble ( xml->GetText() );

		if ( s == "TTmax" )    TTmax    = ToDouble ( xml->GetText() );

		xml = n->GetNextChild ( "DOUBLE" );
	}
}

void Model1v2::init() {

	double WagT[22][21] = {
		{+0.23680e-1,+0.15260e+0,+0.66146e-2,+0.97879e-3,+0.60002e-2,-0.14447e-1,-0.19974e-2,+0.16256e-2,-0.20861e-2,+0.10054e-1,-0.59153e-3,+0.50710e-2,+0.11960e-2,+0.13316e-2,+0.43638e-3,+0.49284e-2,+0.46434e-3,+0.42182e-2,+0.30748e-3,+0.55746e-3,-0.96705e-4},
		{+0.00000e+0,-0.45175e+0,-0.24044e-2,-0.21257e-1,-0.68887e-2,+0.67897e-1,+0.38679e-3,+0.15350e-1,+0.25304e-2,-0.39916e-2,-0.18073e-2,+0.68506e-2,+0.13320e-2,+0.64404e-2,+0.17048e-2,-0.10662e-2,-0.18045e-2,-0.13301e-2,-0.17601e-2,-0.16006e-3,+0.66625e-3},
		{+0.32155e-1,+0.16095e+0,+0.14351e-2,+0.28507e-1,-0.42671e-2,-0.53118e-1,+0.10643e-2,+0.37410e-1,-0.60685e-2,-0.13500e-1,+0.46384e-2,+0.13219e-1,-0.29566e-2,+0.64063e-2,+0.42304e-2,+0.16274e-2,-0.22250e-2,+0.60626e-2,+0.26669e-2,-0.11136e-2,+0.21718e-4},
		{+0.00000e+0,-0.63525e+0,-0.11844e-1,+0.52180e-1,-0.97199e-2,+0.63724e-1,+0.93564e-2,-0.70870e-2,-0.42418e-2,+0.17681e-1,+0.32632e-2,-0.14014e-2,+0.14379e-2,+0.99927e-2,-0.95800e-3,-0.66572e-2,+0.19305e-2,+0.45298e-2,-0.20864e-2,-0.51539e-2,+0.12468e-2},
		{+0.25350e-1,+0.17820e+0,+0.14674e-1,+0.28054e-1,-0.16328e-1,-0.53041e-1,+0.60605e-3,+0.36823e-1,-0.25429e-2,-0.17689e-1,+0.27331e-2,+0.21436e-1,-0.24782e-2,+0.12317e-2,+0.50980e-2,+0.78076e-2,-0.37816e-2,+0.35353e-2,+0.53014e-2,+0.21940e-2,-0.28306e-2},
		{+0.00000e+0,-0.74777e+0,-0.13822e-1,+0.10077e+0,-0.11318e-1,+0.47186e-1,+0.10666e-1,-0.90239e-2,-0.78452e-2,+0.23941e-1,+0.80787e-2,-0.14942e-3,-0.31925e-2,+0.92620e-2,+0.15527e-2,-0.65683e-2,-0.61655e-3,+0.51033e-2,-0.60263e-3,-0.82244e-2,-0.63789e-3},
		{+0.24182e-1,+0.19526e+0,+0.24861e-1,+0.20096e-1,-0.27491e-1,-0.49417e-1,+0.10719e-1,+0.31616e-1,-0.60703e-2,-0.16097e-1,+0.36134e-3,+0.23123e-1,-0.38763e-4,+0.15034e-2,+0.73161e-2,+0.10954e-1,-0.62652e-2,-0.47583e-5,+0.64641e-2,+0.36168e-2,-0.41148e-2},
		{+0.00000e+0,-0.83513e+0,-0.17626e-1,+0.12925e+0,-0.16164e-1,+0.35073e-1,+0.17977e-1,+0.94574e-3,-0.10011e-1,+0.18813e-1,+0.96321e-2,+0.26472e-2,-0.34897e-2,+0.74367e-2,+0.70138e-3,-0.59772e-2,+0.63688e-3,+0.59988e-2,+0.97174e-3,-0.82817e-2,-0.28370e-2},
		{+0.34168e-1,+0.21336e+0,+0.15997e-1,+0.18491e-1,-0.28473e-1,-0.57102e-1,+0.13258e-1,+0.31091e-1,-0.61834e-2,-0.12518e-1,+0.14889e-2,+0.24639e-1,-0.43034e-3,-0.21774e-2,+0.62884e-2,+0.12499e-1,-0.57941e-2,-0.19753e-2,+0.72177e-2,+0.72116e-2,-0.40411e-2},
		{+0.00000e+0,-0.89037e+0,-0.19243e-1,+0.14229e+0,-0.18926e-1,+0.26496e-1,+0.19985e-1,+0.95880e-2,-0.15497e-1,+0.16351e-1,+0.12108e-1,+0.60098e-2,-0.40057e-2,+0.66508e-2,+0.26481e-2,-0.65500e-2,-0.25089e-2,+0.51549e-2,+0.34309e-2,-0.80479e-2,-0.28547e-2},
		{+0.23451e-1,+0.15234e+0,+0.65422e-2,+0.10679e-2,+0.59102e-2,-0.14467e-1,-0.19208e-2,+0.16123e-2,-0.20361e-2,+0.98661e-2,-0.38284e-3,+0.49358e-2,+0.97493e-3,+0.13744e-2,+0.36877e-3,+0.47643e-2,+0.53905e-3,+0.41980e-2,+0.57987e-3,+0.51957e-3,-0.15836e-4},
		{+0.00000e+0,-0.45189e+0,-0.24681e-2,-0.21427e-1,-0.70269e-2,+0.67759e-1,+0.22291e-3,+0.15526e-1,+0.25736e-2,-0.39835e-2,-0.18292e-2,+0.67522e-2,+0.12664e-2,+0.65164e-2,+0.15957e-2,-0.11204e-2,-0.18836e-2,-0.14484e-2,-0.16466e-2,-0.41572e-3,+0.75927e-3},
		{+0.31922e-1,+0.16114e+0,+0.99198e-3,+0.28523e-1,-0.46824e-2,-0.53162e-1,+0.74323e-3,+0.37786e-1,-0.59009e-2,-0.13650e-1,+0.48123e-2,+0.13356e-1,-0.28356e-2,+0.66011e-2,+0.41218e-2,+0.12833e-2,-0.24525e-2,+0.61399e-2,+0.26417e-2,-0.96108e-3,+0.22652e-5},
		{+0.00000e+0,-0.63528e+0,-0.11945e-1,+0.52443e-1,-0.96855e-2,+0.63872e-1,+0.94490e-2,-0.70557e-2,-0.41397e-2,+0.17913e-1,+0.32751e-2,-0.15649e-2,+0.12139e-2,+0.98431e-2,-0.11306e-2,-0.68633e-2,+0.16502e-2,+0.44140e-2,-0.20713e-2,-0.51416e-2,+0.11161e-2},
		{+0.24744e-1,+0.19698e+0,+0.24688e-1,+0.19970e-1,-0.29298e-1,-0.50489e-1,+0.92885e-2,+0.31605e-1,-0.51017e-2,-0.14934e-1,+0.12086e-2,+0.23772e-1,+0.14059e-3,+0.31704e-3,+0.66002e-2,+0.97627e-2,-0.59991e-2,+0.22851e-3,+0.73574e-2,+0.53360e-2,-0.36541e-2},
		{+0.00000e+0,-0.83446e+0,-0.15875e-1,+0.13064e+0,-0.15776e-1,+0.35044e-1,+0.16119e-1,+0.15147e-3,-0.11703e-1,+0.18128e-1,+0.10143e-1,+0.37940e-2,-0.25960e-2,+0.86001e-2,+0.75758e-3,-0.63507e-2,-0.36952e-3,+0.42495e-2,+0.12208e-2,-0.82948e-2,-0.19871e-2},
		{+0.34061e-1,+0.21365e+0,+0.15858e-1,+0.18327e-1,-0.28700e-1,-0.57304e-1,+0.12995e-1,+0.30856e-1,-0.64453e-2,-0.12383e-1,+0.13141e-2,+0.24622e-1,-0.62161e-3,-0.19550e-2,+0.68171e-2,+0.12593e-1,-0.57680e-2,-0.18257e-2,+0.74443e-2,+0.72116e-2,-0.41164e-2},
		{+0.00000e+0,-0.89003e+0,-0.19294e-1,+0.14235e+0,-0.18767e-1,+0.26340e-1,+0.20316e-1,+0.94515e-2,-0.15484e-1,+0.16308e-1,+0.11972e-1,+0.59032e-2,-0.43512e-2,+0.67547e-2,+0.28778e-2,-0.67282e-2,-0.20573e-2,+0.47559e-2,+0.33164e-2,-0.79023e-2,-0.24768e-2},
		{+0.17533e-1,+0.18502e+0,+0.15280e-1,+0.29186e-1,-0.14518e-1,-0.72663e-1,+0.18334e-2,+0.58557e-1,-0.12843e-1,-0.22873e-1,+0.62426e-2,+0.17177e-1,-0.39038e-2,+0.14212e-1,+0.55205e-2,-0.77101e-3,-0.74621e-3,+0.12434e-1,+0.42080e-2,-0.49419e-2,+0.15046e-2},
		{+0.00000e+0,-0.73186e+0,-0.10683e-1,+0.83232e-1,-0.67230e-2,+0.74907e-1,-0.57188e-2,-0.29067e-1,+0.69646e-2,+0.31902e-1,+0.11766e-2,-0.67097e-2,+0.27271e-2,+0.14362e-1,-0.16196e-2,-0.14496e-1,+0.47854e-2,+0.42777e-2,-0.54757e-2,-0.88629e-2,+0.24095e-2},
		{+0.23478e-1,+0.18045e+0,+0.20771e-1,+0.28480e-1,-0.17734e-1,-0.44127e-1,+0.26768e-2,+0.29516e-1,+0.84623e-3,-0.13558e-1,+0.20876e-2,+0.21133e-1,-0.17863e-2,+0.17519e-2,+0.60328e-2,+0.97352e-2,-0.36217e-2,+0.38472e-2,+0.44636e-2,+0.32672e-2,-0.40864e-2},
		{+0.00000e+0,-0.77599e+0,-0.93828e-2,+0.11419e+0,-0.12502e-1,+0.44333e-1,+0.12637e-1,-0.15496e-2,-0.79187e-2,+0.16728e-1,+0.79242e-2,+0.43160e-2,-0.59656e-2,+0.69874e-2,+0.54287e-3,-0.46027e-2,-0.13075e-2,+0.10496e-2,-0.25852e-3,-0.74776e-2,-0.12379e-2}
	};

	double WagQ[26][21] = {
		{+0.32412e-1,+0.22769e-0,-0.27622e-1,+0.41402e-2,+0.12264e-1,-0.13100e-2,+0.24133e-2,-0.72338e-2,-0.47773e-3,+0.16189e-1,-0.53801e-2,+0.14161e-1,-0.46008e-3,+0.25050e-3,-0.19800e-2,+0.57280e-2,-0.40970e-2,+0.75132e-2,+0.50548e-3,-0.57544e-4,+0.12067e-2},
		{+0.00000e-0,-0.62757e-0,+0.35503e-1,-0.41952e-1,-0.87073e-2,+0.11366e-0,-0.63308e-2,+0.31429e-1,-0.51207e-3,-0.12663e-1,-0.20723e-2,+0.10887e-1,-0.24566e-2,+0.10271e-1,+0.15357e-2,-0.68530e-2,-0.11526e-2,-0.20622e-2,-0.13874e-2,+0.14868e-2,+0.25205e-2},
		{+0.47734e-1,+0.24469e-0,-0.34342e-1,+0.48075e-1,+0.14028e-1,-0.47355e-1,-0.13162e-2,+0.20234e-1,-0.26110e-1,+0.10196e-1,+0.20437e-1,-0.13622e-2,-0.12381e-1,+0.24253e-1,-0.59326e-3,-0.68631e-3,+0.55097e-2,+0.66757e-2,-0.51424e-2,+0.29728e-2,+0.68255e-2},
		{+0.00000e-0,-0.88736e-0,+0.23609e-1,+0.54166e-1,-0.31836e-1,+0.12457e-0,+0.14249e-1,-0.27189e-1,-0.12249e-1,+0.26454e-1,+0.46334e-2,+0.13838e-1,+0.68735e-2,+0.58492e-3,-0.13377e-1,-0.11603e-2,+0.13264e-1,-0.42987e-2,-0.76084e-2,+0.77164e-4,-0.37948e-2},
		{+0.24645e-1,+0.26718e-0,+0.16056e-1,+0.65822e-1,-0.22497e-1,-0.78062e-1,+0.24126e-2,+0.61475e-1,-0.16065e-1,-0.33291e-1,+0.12311e-1,+0.31123e-1,-0.12559e-1,+0.13948e-1,+0.88397e-2,+0.50358e-3,-0.79990e-2,+0.13345e-1,+0.11928e-1,-0.13556e-2,-0.70825e-2},
		{+0.00000e-0,-0.11081e+1,+0.15909e-2,+0.16455e-0,-0.20601e-1,+0.85343e-1,+0.87856e-2,-0.31327e-1,-0.96650e-2,+0.43190e-1,+0.12453e-1,+0.95986e-3,-0.79986e-2,+0.15073e-1,+0.24595e-2,-0.16918e-1,+0.51603e-2,+0.11504e-1,-0.47976e-2,-0.14566e-1,+0.23280e-2},
		{+0.38330e-1,+0.30834e-0,+0.13009e-1,+0.37977e-1,-0.42771e-1,-0.62278e-1,+0.11654e-1,+0.46308e-1,-0.20672e-1,-0.25171e-1,+0.93390e-2,+0.40552e-1,-0.36183e-2,+0.36580e-2,+0.79060e-2,+0.83125e-2,-0.12199e-1,+0.53102e-2,+0.13882e-1,+0.35419e-2,-0.10801e-1},
		{+0.00000e-0,-0.12379e+1,-0.76857e-2,+0.21230e-0,-0.34173e-1,+0.51927e-1,+0.25757e-1,-0.13343e-1,-0.16388e-1,+0.36642e-1,+0.17869e-1,+0.54887e-2,-0.72657e-2,+0.52610e-2,-0.32526e-2,-0.10724e-1,+0.58753e-2,+0.90872e-2,-0.33639e-2,-0.15324e-1,+0.15551e-2},
		{+0.55885e-1,+0.33553e-0,-0.55351e-2,+0.29565e-1,-0.39970e-1,-0.59824e-1,+0.87083e-2,+0.38511e-1,-0.21977e-1,-0.16701e-1,+0.74012e-2,+0.41122e-1,+0.14304e-2,+0.24615e-2,+0.54771e-2,+0.99097e-2,-0.11674e-1,+0.28215e-2,+0.13405e-1,+0.61852e-2,-0.10103e-1},
		{+0.00000e-0,-0.13087e+1,-0.13050e-1,+0.22932e-0,-0.39610e-1,+0.35474e-1,+0.30131e-1,+0.22306e-2,-0.19345e-1,+0.31211e-1,+0.24556e-1,+0.82519e-2,-0.10001e-1,+0.29153e-2,-0.14795e-2,-0.93816e-2,+0.31241e-2,+0.82514e-2,+0.68793e-3,-0.14969e-1,-0.15397e-2},
		{+0.32922e-1,+0.22818e-0,-0.27076e-1,+0.40533e-2,+0.11865e-1,-0.16010e-2,+0.16364e-2,-0.75320e-2,-0.10827e-2,+0.16051e-1,-0.54749e-2,+0.14265e-1,-0.85995e-4,+0.40868e-3,-0.12955e-2,+0.59046e-2,-0.38788e-2,+0.73831e-2,+0.41273e-3,-0.54923e-3,+0.10712e-2},
		{+0.00000e-0,-0.62722e-0,+0.35740e-1,-0.41363e-1,-0.81833e-2,+0.11410e-0,-0.58649e-2,+0.31366e-1,-0.64464e-3,-0.13427e-1,-0.25669e-2,+0.10360e-1,-0.27531e-2,+0.10084e-1,+0.18382e-2,-0.67291e-2,-0.70726e-3,-0.17222e-2,-0.10638e-2,+0.17638e-2,+0.27112e-2},
		{+0.47745e-1,+0.24474e-0,-0.33679e-1,+0.48213e-1,+0.14481e-1,-0.47373e-1,-0.23258e-2,+0.19608e-1,-0.26832e-1,+0.10348e-1,+0.20883e-1,-0.78899e-3,-0.12592e-1,+0.24307e-1,-0.43386e-3,-0.40164e-3,+0.60723e-2,+0.69114e-2,-0.47407e-2,+0.22412e-2,+0.62279e-2},
		{+0.00000e-0,-0.88736e-0,+0.23048e-1,+0.54388e-1,-0.30945e-1,+0.12533e-0,+0.15077e-1,-0.27471e-1,-0.12804e-1,+0.25570e-1,+0.43878e-2,+0.13814e-1,+0.69907e-2,+0.69677e-3,-0.14066e-1,-0.14160e-2,+0.13130e-1,-0.35743e-2,-0.68982e-2,+0.60408e-3,-0.35655e-2},
		{+0.38821e-1,+0.30857e-0,+0.12968e-1,+0.37596e-1,-0.43202e-1,-0.62517e-1,+0.11711e-1,+0.46338e-1,-0.20997e-1,-0.25065e-1,+0.94712e-2,+0.40479e-1,-0.34106e-2,+0.34406e-2,+0.72600e-2,+0.84629e-2,-0.12383e-1,+0.55770e-2,+0.14548e-1,+0.31655e-2,-0.10461e-1},
		{+0.00000e-0,-0.12374e+1,-0.69162e-2,+0.21242e-0,-0.34391e-1,+0.51711e-1,+0.25874e-1,-0.12915e-1,-0.16229e-1,+0.36854e-1,+0.17505e-1,+0.53011e-2,-0.66787e-2,+0.51614e-2,-0.27707e-2,-0.11040e-1,+0.50134e-2,+0.87478e-2,-0.33660e-2,-0.14750e-1,+0.19761e-2},
		{+0.55895e-1,+0.33575e-0,-0.54483e-2,+0.30159e-1,-0.39788e-1,-0.59767e-1,+0.85676e-2,+0.38556e-1,-0.22022e-1,-0.16580e-1,+0.76832e-2,+0.40968e-1,+0.14307e-2,+0.19131e-2,+0.56543e-2,+0.99527e-2,-0.11082e-1,+0.27185e-2,+0.13265e-1,+0.63768e-2,-0.10405e-1},
		{+0.00000e-0,-0.13089e+1,-0.12986e-1,+0.22943e-0,-0.39545e-1,+0.35275e-1,+0.29736e-1,+0.19505e-2,-0.19625e-1,+0.31577e-1,+0.24466e-1,+0.80744e-2,-0.10030e-1,+0.27097e-2,-0.14513e-2,-0.94578e-2,+0.27592e-2,+0.82455e-2,+0.74715e-3,-0.15220e-1,-0.14567e-2},
		{+0.36082e-1,+0.27478e-0,-0.60351e-2,+0.70069e-1,-0.79315e-2,-0.10096e-0,-0.11682e-1,+0.80251e-1,-0.17612e-1,-0.25992e-1,+0.14525e-1,+0.13287e-1,-0.13192e-1,+0.31817e-1,+0.11369e-1,-0.91563e-2,-0.12451e-1,+0.20853e-1,+0.13686e-1,-0.66659e-2,-0.74024e-2},
		{+0.00000e-0,-0.10801e+1,+0.40219e-2,+0.12286e-0,-0.28156e-1,+0.12603e-0,+0.16670e-2,-0.56992e-1,-0.39189e-2,+0.48374e-1,+0.14034e-1,-0.66918e-2,-0.11141e-1,+0.20220e-1,+0.53074e-2,-0.20058e-1,+0.66949e-2,+0.57476e-2,-0.75793e-2,-0.93035e-2,+0.84916e-2},
		{+0.33515e-1,+0.28220e-0,+0.10059e-1,+0.48943e-1,-0.24185e-1,-0.64821e-1,+0.70829e-2,+0.62483e-1,-0.88086e-2,-0.32211e-1,+0.74303e-2,+0.29956e-1,-0.11101e-1,+0.98516e-2,+0.99828e-2,+0.99234e-2,-0.13624e-1,+0.37064e-2,+0.81717e-2,+0.37055e-2,-0.58415e-2},
		{+0.00000e-0,-0.11542e+1,+0.11750e-1,+0.18068e-0,-0.30359e-1,+0.79442e-1,+0.21067e-1,-0.17281e-1,-0.20963e-1,+0.25310e-1,+0.14026e-1,+0.79469e-2,-0.10034e-1,+0.11763e-1,+0.19718e-2,-0.13576e-1,-0.29719e-3,+0.60037e-2,+0.16139e-2,-0.13696e-1,-0.30009e-3},
		{+0.00000e-0,-0.80123e-0,-0.83486e-2,+0.12042e-0,-0.14185e-1,+0.43515e-1,+0.13807e-1,+0.42871e-2,-0.62870e-2,+0.12691e-1,+0.77749e-2,+0.63927e-2,-0.49631e-2,+0.44461e-2,+0.49405e-3,-0.32926e-2,-0.22670e-2,-0.62284e-3,+0.62955e-3,-0.68236e-2,-0.19772e-2},
		{+0.00000e-0,-0.11890e+1,+0.22405e-1,+0.18961e-0,-0.35017e-1,+0.81917e-1,+0.28318e-1,-0.53679e-2,-0.22347e-1,+0.15508e-1,+0.10065e-1,+0.14999e-1,-0.56936e-2,+0.80332e-2,+0.23691e-2,-0.82944e-2,-0.51470e-2,+0.10568e-3,+0.59189e-2,-0.11771e-1,-0.29681e-2},
		{+0.00000e-0,-0.82637e-0,-0.44265e-2,+0.12519e-0,-0.15080e-1,+0.45207e-1,+0.12112e-1,+0.84001e-2,-0.36180e-2,+0.85657e-2,+0.50675e-2,+0.85636e-2,-0.54398e-2,+0.40116e-2,+0.15876e-2,-0.25855e-2,-0.39788e-2,-0.77822e-3,+0.14790e-2,-0.60376e-2,-0.26079e-2},
		{+0.00000e-0,-0.12213e+1,+0.27934e-1,+0.18974e-0,-0.35626e-1,+0.86348e-1,+0.32399e-1,+0.36867e-2,-0.17751e-1,+0.79699e-2,+0.25833e-2,+0.14951e-1,-0.27002e-2,+0.81931e-2,+0.36117e-2,-0.49921e-2,-0.64297e-2,-0.30697e-2,+0.62224e-2,-0.10872e-1,-0.28107e-2}
	};


	l     = L / 2; 			//Half-length (m)
	d     = ( dfor + daft ) / 2; 		//Draft (m)
	m     = rho * Cb * LWL * d * BOA; 	//Mass (kg)
	Izz   = m * pow ( Rg,2 ); 		//Moment of inertia
	AP    = M_PI * pow ( DP/2,2 ); 		//Area of propellor (m^2)
	HR    = AR/ChordR;			//Height of rudder (m)
	eta   = DP/HR;			//Ratio diameter of propeller - height of rudder (-)

	prtT = 0;

	for ( int i=0; i<11; i++ ) prtT += prT[i];

	for ( int j=0; j<21; j++ ) {
		AKT[j] = 0;
		BKT[j] = 0;
	}

	for ( int i=0; i<11; i++ ) {
		for ( int j=0; j<21; j++ ) {
			AKT[j] += WagT[2*i  ][j]*prT[i]/prtT;
			BKT[j] += WagT[2*i+1][j]*prT[i]/prtT;
		}
	}

	prtQ = 0;

	for ( int i=0; i<11; i++ ) prQ[i] = 1.;

	for ( int i=0; i<11; i++ ) prtQ += prQ[i];

	for ( int j=0; j<21; j++ ) {
		AKQ[j] = 0;
		BKQ[j] = 0;
	}

	for ( int i=0; i<11; i++ ) {
		for ( int j=0; j<21; j++ ) {
			AKQ[j] += WagQ[2*i  ][j]*prQ[i]/prtQ;
			BKQ[j] += WagQ[2*i+1][j]*prQ[i]/prtQ;
		}
	}

	xP = ndxP * L;
	yP = ndyP * L;

	//Ideal fluid effects
	Xup = ndXup * m;
	Xrr = ndXrr * m * L;
	Yvp = ndYvp * m;
	Nrp = ndNrp * m * L * L;
	Xvr = ndXvr * m;
	Xvv = ndXvv * m / L;
	Yrp = ndYrp * m * L;
	Nvp = ndNvp * m * L;

	//Hull resistance effects
	XHRu   = ndXHRu   * m * sqrt ( g / L );
	XHRuau = ndXHRuau * m / L;
	XHRuuu = ndXHRuuu * m / L / sqrt ( L * g );

	YPTp = ndYPTp;
	YPTm = ndYPTm;
	NPTp = ndNPTp * L;
	NPTm = ndNPTm * L;

	xR = ndxR * L;	//Longitudinal position of rudder forward of midship (m)

	m11 = m - Xup;
	m22 = m - Yvp;
	m23 = m*xG - Yrp;
	m32 = m*xG - Nvp;
	m33 = Izz - Nrp;
	det = m22*m33 - m32*m23;
}

void Model1v2::get_bounds ( double *low, double *upp ) {
	low[0]  =  nan ( "" );
	low[1]  =  nan ( "" );
	low[2]  =  nan ( "" );
	low[3]  =  nan ( "" );
	low[4]  =  nan ( "" );
	low[5]  =  nan ( "" );
	low[6]  =  nan ( "" );
	low[7]  =  nan ( "" );
	low[8]  =  nan ( "" ); //0.1;//-1.5;
	low[9]  =  nan ( "" ); //0.1;//-1.5;
//     low[10] =  -35.*M_PI/180;
//     low[11] =  -35.*M_PI/180;
//     low[12] =  -1.;
//     low[13] =  -1.;


	upp[0]  = nan ( "" );
	upp[1]  = nan ( "" );
	upp[2]  = nan ( "" );
	upp[3]  = nan ( "" );
	upp[4]  = nan ( "" );
	upp[5]  = nan ( "" );
	upp[6]  = nan ( "" );
	upp[7]  = nan ( "" );
	upp[8]  = nan ( "" ); //+1.5;
	upp[9]  = nan ( "" ); //+1.5;
//     upp[10] = +35.*M_PI/180;
//     upp[11] = +35.*M_PI/180;
//     upp[12] = +1.;
//     upp[13] = +1.;

//     low[0] = -35.*M_PI/180;
//     low[1] = -35.*M_PI/180;
// //     low[2] = -98.8/60;
// //     low[3] = -98.8/60;
//     low[2] = -1;
//     low[3] = -1;
// //     low[4] = -1;
//
//     upp[0] = +35.*M_PI/180;
//     upp[1] = +35.*M_PI/180;
// //     upp[2] = +98.8/60;
// //     upp[3] = +98.8/60;
//     upp[2] = +1;
//     upp[3] = +1;
// //     upp[4] = +1;
}

void Model1v2::get_diff_bounds ( double *low, double *upp ) {
//     low[0] = -3.*M_PI/180;
//     low[1] = -3.*M_PI/180;
//     low[2] = -0.01;
//     low[3] = -0.01;
// //     low[4] = -1;
//
//     upp[0] = +3.*M_PI/180;
//     upp[1] = +3.*M_PI/180;
//     upp[2] = +0.01;
//     upp[3] = +0.01;
// //     upp[4] = +1;
}

int Model1v2::get_nbrStates() {
//     return 8;
	return 10;
}

int Model1v2::get_nbrControls() {
	return 0;
}













