#include "functions.h"
#include "../global.h"

void convert(string datatype, int *v) {
	if(datatype == "uInt8" ) *v = (uint8_t )(*v);
	if(datatype == "Int8"  ) *v = (int8_t  )(*v);
	if(datatype == "uInt16") *v = (uint16_t)(*v);
	if(datatype == "Int16" ) *v = (int16_t )(*v);
	if(datatype == "uInt32") *v = (uint32_t)(*v);
	if(datatype == "Int32" ) *v = (int32_t )(*v);
}

double min2(double a, double b) {
	if(a<=b) {
		return a;
	} else {
		return b;
	}
}

// sin in degrees
double sind(double x){
	return sin(fmod((x),360) * M_PI / 180);
}

//cos in degrees
double cosd(double x){
	return cos(fmod((x),360) * M_PI / 180);
}

void LLAtoNED(double lat, double lon, double h, double lat0, double lon0, double h0,
			  double *x, double *y, double *z) {
	
	// calculation taken from:
	//
	// Unmanned Rotorcraft Systems
	// Cai, Guowei, Chen, Ben M., Lee, Tong Heng
	// Advances in Industrial Control, 2011, Springer-Verlag London
	// chap. 2, pp. 23-32

	// transformation from LLA to NED
	double R_ea=6378137.0; //[m]
    double f=1/298.257223563;
    double R_eb=R_ea*(1.-f); //[m]
    double e=sqrt(R_ea*R_ea-R_eb*R_eb)/R_ea;
    
    double N_e=R_ea/sqrt(1.-e*e*sind(lat)*sind(lat));
    double N_e0=R_ea/sqrt(1.-e*e*sind(lat0)*sind(lat0));
    
    double P_e_1=(N_e+h)*cosd(lat)*cosd(lon) - (N_e0+h0)*cosd(lat0)*cosd(lon0);
	double P_e_2=(N_e+h)*cosd(lat)*sind(lon) - (N_e0+h0)*cosd(lat0)*sind(lon0);
	double P_e_3=(N_e*(1.-e*e)+h)*sind(lat)  - (N_e0*(1-e*e)+h0)*sind(lat0);
	
	// matrix-multiplication
	*x = -sind(lat0)*cosd(lon0)*P_e_1 + -sind(lat0)*sind(lon0)*P_e_2 +  cosd(lat0)*P_e_3;
	*y = -sind(lon0)           *P_e_1 +  cosd(lon0)           *P_e_2 +  0.        *P_e_3;
	*z = -cosd(lat0)*cosd(lon0)*P_e_1 + -cosd(lat0)*sind(lon0)*P_e_2 + -sind(lat0)*P_e_3;

}

void NEDtoLLA(double x, double y, double z, double lat0, double lon0, double h0, 
			  double *lat, double *lon, double *h) {
	
	// calculation taken from:
	//
	// Unmanned Rotorcraft Systems
	// Cai, Guowei, Chen, Ben M., Lee, Tong Heng
	// Advances in Industrial Control, 2011, Springer-Verlag London
	// chap. 2, pp. 23-32

	// transformation from NED to LLA
	double R_ea=6378137.0; //[m]
    double f=1/298.257223563;
    double R_eb=R_ea*(1.-f); //[m]
    double e=sqrt(R_ea*R_ea-R_eb*R_eb)/R_ea;
	double ep=sqrt(R_ea*R_ea-R_eb*R_eb)/R_eb;
    
    double N_e0=R_ea/sqrt(1.-e*e*sind(lat0)*sind(lat0));
	
    double P_e0_1=(N_e0+h0)*cosd(lat0)*cosd(lon0);
	double P_e0_2=(N_e0+h0)*cosd(lat0)*sind(lon0);
	double P_e0_3=(N_e0*(1.-e*e)+h0)*sind(lat0);
	
	double x_e = -sind(lat0)*cosd(lon0)*x-sind(lon0)*y-cosd(lat0)*cosd(lon0)*z + P_e0_1;
	double y_e = -sind(lat0)*sind(lon0)*x+cosd(lon0)*y-cosd(lat0)*sind(lon0)*z + P_e0_2;
	double z_e = cosd(lat0)*x+0.*y-sind(lat0)*z + P_e0_3;
	
	*lon = atan2(y_e, x_e)*180./M_PI;
	
	double p=sqrt(x_e*x_e+y_e*y_e);
    double theta=atan2(z_e*R_ea,p*R_eb)*180./M_PI;
	
	*lat = atan2(z_e+ep*ep*R_eb*sind(theta)*sind(theta)*sind(theta),p-e*e*R_ea*cosd(theta)*cosd(theta)*cosd(theta))*180/M_PI;
    
    double N_e=R_ea/sqrt(1-e*e*sind(*lat)*sind(*lat));
    
    *h= p/cosd(*lat)-N_e;
    
}

unsigned nChoosek ( unsigned n, unsigned k ) {
	if ( k > n ) return 0;

	if ( k * 2 > n ) k = n-k;

	if ( k == 0 ) return 1;

	int result = n;

	for ( int i = 2; i <= k; ++i ) {
		result *= ( n-i+1 );
		result /= i;
	}

	return result;
}


void cubicSpline(int i, int k, double dt, double dtfull, bool extrapolate, double *x0, double *y0, double *psi) {
	double x1;
	double x2;
	double xp1;
	double xp2;
	double a0;
	double a1;
	double a2;
	double a3;
	
	if(!extrapolate) {
		if(dt !=0) {
			x1  = ph0->x(i  , ph0->NS[k]);
			x2  = ph0->x(i+1, ph0->NS[k]);
			xp1 = ph0->scale[k][3]*ph0->x(i  , ph0->NS[k]+3) * cos ( ph0->scale[k][2]*ph0->x(i  , ph0->NS[k]+2) ) - ph0->scale[k][4]*ph0->x(i  , ph0->NS[k]+4) * sin ( ph0->scale[k][2]*ph0->x(i  , ph0->NS[k]+2) );
			xp2 = ph0->scale[k][3]*ph0->x(i+1, ph0->NS[k]+3) * cos ( ph0->scale[k][2]*ph0->x(i+1, ph0->NS[k]+2) ) - ph0->scale[k][4]*ph0->x(i+1, ph0->NS[k]+4) * sin ( ph0->scale[k][2]*ph0->x(i+1, ph0->NS[k]+2) );

			xp1 *= dtfull/ph0->scale[k][0];
			xp2 *= dtfull/ph0->scale[k][0];
			
			a0  = xp1+xp2-2*(x2-x1);
			a1  = 3*(x2-x1)-2*xp1-xp2;
			a2  = xp1;
			a3  = x1;
			
			*x0 = a0*pow(dt,3) + a1*pow(dt,2) + a2*dt + a3;
			
			x1  = ph0->x(i  , ph0->NS[k]+1);
			x2  = ph0->x(i+1, ph0->NS[k]+1);
			xp1 = ph0->scale[k][3]*ph0->x(i  , ph0->NS[k]+3) * sin ( ph0->scale[k][2]*ph0->x(i  , ph0->NS[k]+2) ) + ph0->scale[k][4]*ph0->x(i  , ph0->NS[k]+4) * cos ( ph0->scale[k][2]*ph0->x(i  , ph0->NS[k]+2) );
			xp2 = ph0->scale[k][3]*ph0->x(i+1, ph0->NS[k]+3) * sin ( ph0->scale[k][2]*ph0->x(i+1, ph0->NS[k]+2) ) + ph0->scale[k][4]*ph0->x(i+1, ph0->NS[k]+4) * cos ( ph0->scale[k][2]*ph0->x(i+1, ph0->NS[k]+2) );

			xp1 *= dtfull/ph0->scale[k][1];
			xp2 *= dtfull/ph0->scale[k][1];
			
			a0  = xp1+xp2-2*(x2-x1);
			a1  = 3*(x2-x1)-2*xp1-xp2;
			a2  = xp1;
			a3  = x1;
			
			*y0 = a0*pow(dt,3) + a1*pow(dt,2) + a2*dt + a3;
			
			x1  = ph0->x(i  , ph0->NS[k]+2);
			x2  = ph0->x(i+1, ph0->NS[k]+2);
			xp1 = ph0->scale[k][5]*ph0->x(i  , ph0->NS[k]+5);
			xp2 = ph0->scale[k][5]*ph0->x(i+1, ph0->NS[k]+5);

			xp1 *= dtfull/ph0->scale[k][2];
			xp2 *= dtfull/ph0->scale[k][2];
			
			a0  = xp1+xp2-2*(x2-x1);
			a1  = 3*(x2-x1)-2*xp1-xp2;
			a2  = xp1;
			a3  = x1;
			
			*psi = a0*pow(dt,3) + a1*pow(dt,2) + a2*dt + a3;
			
		} else {
			*x0  = ph0->x(i  , ph0->NS[k]  );
			*y0  = ph0->x(i  , ph0->NS[k]+1);
			*psi = ph0->x(i  , ph0->NS[k]+2);
		}
	} else {
		xp1  = ph0->scale[k][3]*ph0->x(i , ph0->NS[k]+3) * cos ( ph0->scale[k][2]*ph0->x(i , ph0->NS[k]+2) ) - ph0->scale[k][4]*ph0->x(i , ph0->NS[k]+4) * sin ( ph0->scale[k][2]*ph0->x(i , ph0->NS[k]+2) );

		xp1 *= dtfull/ph0->scale[k][0];
		
		*x0  = ph0->x(i , ph0->NS[k]) + dt*xp1;

		xp1  = ph0->scale[k][3]*ph0->x(i , ph0->NS[k]+3) * sin ( ph0->scale[k][2]*ph0->x(i , ph0->NS[k]+2) ) + ph0->scale[k][4]*ph0->x(i , ph0->NS[k]+4) * cos ( ph0->scale[k][2]*ph0->x(i , ph0->NS[k]+2) );

		xp1 *= dtfull/ph0->scale[k][1];
		
		*y0  = ph0->x(i , ph0->NS[k]+1) + dt*xp1;

		*psi = ph0->x(i , ph0->NS[k]+2);
	}
}



void cubicSpline(int i, int k, double dt, double dtfull, bool extrapolate, double *x0, double *y0, double *psi, double *delta, double *n) {
	
	double x1;
	double x2;
	double xp1;
	double xp2;
	double a0;
	double a1;
	double a2;
	double a3;
	
	if(!extrapolate) {
		if(dt !=0) {
			x1  = ph0->scale[k][0]*ph0->x(i  , ph0->NS[k]+0);
			x2  = ph0->scale[k][0]*ph0->x(i+1, ph0->NS[k]+0);
			xp1 = (ph0->scale[k][3]*ph0->x(i  , ph0->NS[k]+3) * cos ( ph0->scale[k][2]*ph0->x(i  , ph0->NS[k]+2) ) - ph0->scale[k][4]*ph0->x(i  , ph0->NS[k]+4) * sin ( ph0->scale[k][2]*ph0->x(i  , ph0->NS[k]+2) ));
			xp2 = (ph0->scale[k][3]*ph0->x(i+1, ph0->NS[k]+3) * cos ( ph0->scale[k][2]*ph0->x(i+1, ph0->NS[k]+2) ) - ph0->scale[k][4]*ph0->x(i+1, ph0->NS[k]+4) * sin ( ph0->scale[k][2]*ph0->x(i+1, ph0->NS[k]+2) ));
			
			xp1 *= dtfull;
			xp2 *= dtfull;
			
			a0  = xp1+xp2-2*(x2-x1);
			a1  = 3*(x2-x1)-2*xp1-xp2;
			a2  = xp1;
			a3  = x1;
			
			*x0 = a0*pow(dt,3) + a1*pow(dt,2) + a2*dt + a3;
			
			x1  = ph0->scale[k][1]*ph0->x(i  , ph0->NS[k]+1);
			x2  = ph0->scale[k][1]*ph0->x(i+1, ph0->NS[k]+1);
			xp1 = (ph0->scale[k][3]*ph0->x(i  , ph0->NS[k]+3) * sin ( ph0->scale[k][2]*ph0->x(i  , ph0->NS[k]+2) ) + ph0->scale[k][4]*ph0->x(i  , ph0->NS[k]+4) * cos ( ph0->scale[k][2]*ph0->x(i  , ph0->NS[k]+2) ));
			xp2 = (ph0->scale[k][3]*ph0->x(i+1, ph0->NS[k]+3) * sin ( ph0->scale[k][2]*ph0->x(i+1, ph0->NS[k]+2) ) + ph0->scale[k][4]*ph0->x(i+1, ph0->NS[k]+4) * cos ( ph0->scale[k][2]*ph0->x(i+1, ph0->NS[k]+2) ));
			
			xp1 *= dtfull/*/ph0->scale[k][1]*/;
			xp2 *= dtfull/*/ph0->scale[k][1]*/;
			
			a0  = xp1+xp2-2*(x2-x1);
			a1  = 3*(x2-x1)-2*xp1-xp2;
			a2  = xp1;
			a3  = x1;
			
			*y0 = a0*pow(dt,3) + a1*pow(dt,2) + a2*dt + a3;
			
			x1  = ph0->scale[k][2]*ph0->x(i  , ph0->NS[k]+2);
			x2  = ph0->scale[k][2]*ph0->x(i+1, ph0->NS[k]+2);
			xp1 = ph0->scale[k][5]*ph0->x(i  , ph0->NS[k]+5);
			xp2 = ph0->scale[k][5]*ph0->x(i+1, ph0->NS[k]+5);
			
			xp1 *= dtfull;
			xp2 *= dtfull;
			
			a0  = xp1+xp2-2*(x2-x1);
			a1  = 3*(x2-x1)-2*xp1-xp2;
			a2  = xp1;
			a3  = x1;
			
			*psi = a0*pow(dt,3) + a1*pow(dt,2) + a2*dt + a3;
			
			x1  = ph0->scale[k][6]*ph0->x(i  , ph0->NS[k]+6);
			x2  = ph0->scale[k][6]*ph0->x(i+1, ph0->NS[k]+6);
			xp1 = ph0->scale_ctrl[k][0]*ph0->u(i  , ph0->NC[k]+0)*dtfull;
			xp2 = ph0->scale_ctrl[k][0]*ph0->u(i+1, ph0->NC[k]+0)*dtfull;
			a0  = xp1+xp2-2*(x2-x1);
			a1  = 3*(x2-x1)-2*xp1-xp2;
			a2  = xp1;
			a3  = x1;
			
			*delta = a0*pow(dt,3) + a1*pow(dt,2) + a2*dt + a3;
			
			x1  = ph0->scale[k][7]*ph0->x(i  , ph0->NS[k]+7);
			x2  = ph0->scale[k][7]*ph0->x(i+1, ph0->NS[k]+7);
			xp1 = ph0->scale_ctrl[k][1]*ph0->u(i  , ph0->NC[k]+1)*dtfull;
			xp2 = ph0->scale_ctrl[k][1]*ph0->u(i+1, ph0->NC[k]+1)*dtfull;
			a0  = xp1+xp2-2*(x2-x1);
			a1  = 3*(x2-x1)-2*xp1-xp2;
			a2  = xp1;
			a3  = x1;
			
			*n = a0*pow(dt,3) + a1*pow(dt,2) + a2*dt + a3;
		} else {
			*x0    = ph0->scale[k][0]*ph0->x(i  , ph0->NS[k]+0);
			*y0    = ph0->scale[k][1]*ph0->x(i  , ph0->NS[k]+1);
			*psi   = ph0->scale[k][2]*ph0->x(i  , ph0->NS[k]+2);
			*delta = ph0->scale[k][6]*ph0->x(i  , ph0->NS[k]+6);
			*n     = ph0->scale[k][7]*ph0->x(i  , ph0->NS[k]+7);
		}
	} else {
		xp1  = ph0->scale[k][3]*ph0->x(i , ph0->NS[k]+3) * cos ( ph0->scale[k][2]*ph0->x(i , ph0->NS[k]+2) ) - ph0->scale[k][4]*ph0->x(i , ph0->NS[k]+4) * sin ( ph0->scale[k][2]*ph0->x(i , ph0->NS[k]+2) );
		*x0  = ph0->scale[k][0]*ph0->x(i , ph0->NS[k]+0) + dt*xp1;

		xp1  = ph0->scale[k][3]*ph0->x(i , ph0->NS[k]+3) * sin ( ph0->scale[k][2]*ph0->x(i , ph0->NS[k]+2) ) + ph0->scale[k][4]*ph0->x(i , ph0->NS[k]+4) * cos ( ph0->scale[k][2]*ph0->x(i , ph0->NS[k]+2) );
		*y0  = ph0->scale[k][1]*ph0->x(i , ph0->NS[k]+1) + dt*xp1;

		*psi   = ph0->scale[k][2]*ph0->x(i , ph0->NS[k]+2);
		*delta = ph0->scale[k][6]*ph0->x(i , ph0->NS[k]+6);
		*n     = ph0->scale[k][7]*ph0->x(i , ph0->NS[k]+7);
	}
}


void cubicSpline_structure(tw::DiffStructure &s, int i, int k, int idx, CubicSplineBool disc) {
	if(disc == NORMAL) {
		s ( idx, ph0->x_index ( i  ,ph0->NS[k]   ) );
		s ( idx, ph0->x_index ( i  ,ph0->NS[k]+1 ) );
		s ( idx, ph0->x_index ( i  ,ph0->NS[k]+2 ) );
		s ( idx, ph0->x_index ( i  ,ph0->NS[k]+3 ) );
		s ( idx, ph0->x_index ( i  ,ph0->NS[k]+4 ) );
		s ( idx, ph0->x_index ( i  ,ph0->NS[k]+5 ) );
		s ( idx, ph0->x_index ( i+1,ph0->NS[k]   ) );
		s ( idx, ph0->x_index ( i+1,ph0->NS[k]+1 ) );
		s ( idx, ph0->x_index ( i+1,ph0->NS[k]+2 ) );
		s ( idx, ph0->x_index ( i+1,ph0->NS[k]+3 ) );
		s ( idx, ph0->x_index ( i+1,ph0->NS[k]+4 ) );
		s ( idx, ph0->x_index ( i+1,ph0->NS[k]+5 ) );
	}
	if(disc == SINGLEPOINT) {
		s ( idx, ph0->x_index ( i     ,ph0->NS[k]   ) );
		s ( idx, ph0->x_index ( i     ,ph0->NS[k]+1 ) );
		s ( idx, ph0->x_index ( i     ,ph0->NS[k]+2 ) );
	}
	if(disc == EXTRAPOLATE) {
		s ( idx, ph0->x_index ( i  ,ph0->NS[k]   ) );
		s ( idx, ph0->x_index ( i  ,ph0->NS[k]+1 ) );
		s ( idx, ph0->x_index ( i  ,ph0->NS[k]+2 ) );
		s ( idx, ph0->x_index ( i  ,ph0->NS[k]+3 ) );
		s ( idx, ph0->x_index ( i  ,ph0->NS[k]+4 ) );
	}
}


double t_max() {
	double max = 0;
	if(ph0->ENDTIME == OPEN) {
		max = ph0->p(0)*ph0->scale_t[0];
		for(int i=0; i<ph0->nbrCShips; i++) {
			if(ph0->p(i)*ph0->scale_t[i] > max*ph0->scale_t[i]) max = ph0->p(i)*ph0->scale_t[i];
		}
	} 
	if(ph0->ENDTIME == OPEN_WITH_FIXED_RATIOS) {
		max = ph0->p(0)*ph0->scale_t[0]*ph0->ratioToFirst[0];
		for(int i=0; i<ph0->nbrCShips; i++) {
			if(ph0->p(0)*ph0->scale_t[0]*ph0->ratioToFirst[i] > max*ph0->scale_t[0]) max = ph0->p(0)*ph0->scale_t[0]*ph0->ratioToFirst[i];
		}
	} 
	if(ph0->ENDTIME == COMMON_AND_FIXED) {
		max = ph0->tEnd[0]*ph0->scale_t[0];
	}
	return max;
}

void SetFloatTime(tw::Viewer *viewer, double t_scale) {
	//double max = t_max();
	double max = 4000.0;
	viewer->setFloatTime ( max/t_scale, 0 );
}

BaseModel* get_ship_model(string ID) {
	if      (ID.compare("1") == 0) return new Model1v5MeckPomm ( "ParamXML/s1v1_m1v5_p1MeckPomm.xml" );
	else if (ID.compare("2") == 0) return new Model1v5MeckPomm ( "ParamXML/s1v1_m1v5_p1MeckPomm.xml" );
	else return nullptr;
}

BaseModel* get_buggy_model(string ID) {
	if      (ID.compare("1") == 0) return new Buggy_v1 ( "ParamXML/BuggyParams.xml" );
	else if (ID.compare("2") == 0) return new Buggy_v1 ( "ParamXML/BuggyParams.xml" );
	else return nullptr;
}

double scaleplot (int &len, int &ndgl, int &nsteuer, double *t,double *x, double *u, int &j, int &index) {
	int idx = 0;
	for ( int k=0; k<ph0->nbrCShips; k++ ) {
		for ( int i=0; i<ph0->NmodelS[k]; i++ ) {
			if (index == idx ) return ph0->scale[k][i]*ph0->x__(j,ph0->NS[k]+i);
			idx++;
		}
		for ( int i=0; i<ph0->NmodelC[k]; i++ ) {
			if (index == idx ) return ph0->scale[k][ph0->NmodelS[k]+i]*ph0->x__(j,ph0->NS[k]+ph0->NmodelS[k]+i);
			idx++;
		}
		if (index == idx ) return ph0->scale[k][ph0->NmodelS[k]+ph0->NmodelC[k]]*ph0->x__(j,ph0->NS[k]+ph0->NmodelS[k]+ph0->NmodelC[k]);
		idx++;
	}
// 	for ( int i=0; i<ph0->NC[ph0->nbrCShips]; i++ ) {
// 		if (index == idx ) return ph0->u__(j,i);
// 		idx++;
// 	}
	
	for ( int k=0; k<ph0->nbrCShips; k++ ) {
		for ( int i=0; i<ph0->NmodelC[k]; i++ ) {
			if (index == idx ) return ph0->scale_ctrl[k][i]*ph0->u__(j,ph0->NC[k]+i);
			idx++;
		}
	}		
	return 0;
}









// double Clamp( double n, double min, double max ) {
// 	if ( n < min ) return min;
// 
// 	if ( n > max ) return max;
// 
// 	return n;
// }
// 
// double sqDinstanceSegmentSegment( Vec2d &p1, Vec2d &q1, Vec2d &p2, Vec2d &q2 ) {
// 	Vec2d d1 = q1 - p1; // Direction vector of segment S1
// 	Vec2d d2 = q2 - p2; // Direction vector of segment S2
// 	Vec2d r = p1 - p2;
// 	double a = d1.dot( d1 );  // Squared length of segment S1, always nonnegative
// 	double e = d2.dot( d2 );  // Squared length of segment S2, always nonnegative
// 	double f = d2.dot( r );
// 	double s=0.;
// 	double t=0.;
// 
// 	// Check if either or both segments degenerate into points
// 	if ( a <= EPSILON && e <= EPSILON ) {
// 		// Both segments degenerate into points
// 		return ( p1-p2 ).dot( p1-p2 );
// 	}
// 
// 	if ( a <= EPSILON ) {
// 		// First segment degenerates into a point
// 		s = 0.0;
// 		t = f / e; // s = 0 => t = (b*s + f) / e = f / e
// 		t = Clamp( t, 0.0, 1.0 );
// 	} else {
// 		double c = d1.dot( r );
// 
// 		if ( e <= EPSILON ) {
// 			// Second segment degenerates into a point
// 			t = 0.0;
// 			s = Clamp( -c / a, 0.0, 1.0 );  // t = 0 => s = (b*t - c) / a = -c / a
// 		} else {
// 			// The general nondegenerate case starts here
// 			double b = d1.dot( d2 );
// 			double denom = a*e-b*b; // Always nonnegative
// 
// 			// If segments not parallel, compute closest point on L1 to L2 and
// 			// clamp to segment S1. Else pick arbitrary s (here 0)
// 			if ( denom != 0.0 ) {
// 				s = Clamp( ( b*f - c*e ) / denom, 0.0, 1.0 );
// 			} else s = 0.0;
// 
// 			// Compute point on L2 closest to S1(s) using
// 			// t = Dot((P1 + D1*s) - P2,D2) / Dot(D2,D2) = (b*s + f) / e
// 			double tnom = b*s + f;
// 
// 			if ( tnom < 0.0 ) {
// 				t = 0.0;
// 				s = Clamp( -c / a, 0.0, 1.0 );
// 			} else if ( tnom > e ) {
// 				t = 1.0;
// 				s = Clamp( ( b - c ) / a, 0.0, 1.0 );
// 			} else {
// 				t = tnom / e;
// 			}
// 		}
// 	}
// 
// 	return ( ( p1 + d1 * s )- ( p2 + d2 * t ) ).dot( ( p1 + d1 * s )- ( p2 + d2 * t ) );
// }
// 
// // double Area( Poly poly ) {
// // 	double area=0.;
// // 	if ( poly.nbrVertices>2 ) {
// // 		for ( unsigned int i=0; i<poly.nbrVertices-1; i++ ) {
// // 			area=area+( ( poly.vertices_ned[i].x * poly.vertices_ned[i+1].y ) - ( poly.vertices_ned[i].y * poly.vertices_ned[i+1].x ) ); 
// // 		}
// // 		area=area+( ( poly.vertices_ned[poly.nbrVertices-1].x * poly.vertices_ned[0].y ) - ( poly.vertices_ned[poly.nbrVertices-1].y * poly.vertices_ned[0].x ) );
// // 		area=fabs( area )/2.;
// // 	}
// // 	return area;
// // }
// 
// double areaOfPoly( std::vector<Vec2d>& segs ) {
// 	double A=0.;
// 
// 	if ( segs.size()>2 ) {
// 		for ( unsigned int i=0; i<segs.size()-1; i++ ) {
// 			A=A+( ( segs[i].x * segs[i+1].y ) - ( segs[i].y * segs[i+1].x ) ); 
// 		}
// 		A=A+( ( segs[segs.size()-1].x * segs[0].y ) - ( segs[segs.size()-1].y * segs[0].x ) );
// 		A=fabs( A )/2.;
// 	}
// 	return A;
// }
// 
// int line_sect( Vec2d& x0, Vec2d& x1, Vec2d& y0, Vec2d& y1, Vec2d& res ) {
// 	Vec2d dx=x1-x0;
// 	Vec2d dy=y1-y0;
// 	Vec2d d=x0-y0;
// 	/* x0 + a dx = y0 + b dy ->
// 	   x0 X dx = y0 X dx + b dy X dx ->
// 	   b = (x0 - y0) X dx / (dy X dx) */
// 	double dyx = dy.cross( dx );
// 
// 	if ( !dyx ) return 0;
// 
// 	dyx = d.cross( dx ) / dyx;
// 
// 	if ( dyx <= 0 || dyx >= 1 ) return 0;
// 
// 	res.x = y0.x + dyx * dy.x;
// 	res.y = y0.y + dyx * dy.y;
// 	return 1;
// }
// 
// /* tells if vec c lies on the left side of directed edge a->b
//  * 1 if left, -1 if right, 0 if colinear
//  */
// int left_of( Vec2d& a, Vec2d& b, Vec2d& c ) {
// 	double out = ( b-a ).cross( c-b );
// 	return out < 0 ? -1 : out > 0;
// }
// 
// int poly_winding( std::vector<Vec2d>& p ) {
// 	return left_of( p[0], p[1], p[2] );
// }
// 
// void poly_edge_clip( std::vector<Vec2d>& sub, Vec2d& x0, Vec2d& x1, int left, std::vector<Vec2d>& res ) {
// 	int i, side0, side1;
// 	Vec2d tmp;
// 	Vec2d v0 = sub.back();
// 	res.resize( 0 );
// 
// 	side0 = left_of( x0, x1, v0 );
// 
// 	if ( side0 != -left ) res.push_back( v0 );
// 
// 	for ( i = 0; i < sub.size(); i++ ) {
// 		side1 = left_of( x0, x1, sub[i] );
// 
// 		if ( side0 + side1 == 0 && side0 )
// 
// 			/* last point and current straddle the edge */
// 			if ( line_sect( x0, x1, v0, sub[i], tmp ) )
// 				res.push_back( tmp );
// 
// 		if ( i == sub.size()-1 ) break;
// 
// 		if ( side1 != -left ) res.push_back( sub[i] );
// 
// 		v0 = sub[i];
// 		side0 = side1;
// 	}
// }
// 
// void poly_clip( std::vector<Vec2d>& sub, std::vector<Vec2d>& clip, std::vector<Vec2d>& p2) {
// 	std::vector<Vec2d> p1;
// 	p2.resize(0);
// 	std::vector<Vec2d> tmp;
// 
// 	int dir = poly_winding( clip );
// 	poly_edge_clip( sub, clip.back(), clip[0], dir, p2 );
// 
// 	for ( int i = 0; i < clip.size()-1; i++ ) {
// 		tmp = p2;
// 		p2 = p1;
// 		p1 = tmp;
// 
// 		if ( p1.size() == 0 ) {
// 			p2.resize( 0 );
// 			break;
// 		}
// 
// 		poly_edge_clip( p1, clip[i], clip[i+1], dir, p2 );
// 	}
// // 	for(int i=0; i<p2.size(); i++) cout << p2[i] << endl;
// // 	cout << endl;
// }