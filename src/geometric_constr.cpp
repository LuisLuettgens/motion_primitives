#include "geometric_constr.h"

using namespace std;


GeometricConstr::GeometricConstr()
	:
	nbrCShips( 0 ),
	lat_max(-std::numeric_limits<double>::infinity()),
	lat_min( std::numeric_limits<double>::infinity()),
	lon_max(-std::numeric_limits<double>::infinity()),
	lon_min( std::numeric_limits<double>::infinity())
{}

GeometricConstr::GeometricConstr(int nbrCShips)
	:
	nbrCShips( nbrCShips ),
	Cships   ( std::vector<ship> (nbrCShips) ),
	n_harb_polygons( 0 ),
	n_obst_polygons( 0 ),
	n_spee_polygons( 0 ),
	n_guid_polygons( 0 ),
	lat_max(-std::numeric_limits<double>::infinity()),
	lat_min( std::numeric_limits<double>::infinity()),
	lon_max(-std::numeric_limits<double>::infinity()),
	lon_min( std::numeric_limits<double>::infinity())
{
	for(int i=0; i<nbrCShips; i++) {
		Cships[i].polygon.n_segments = n_ship_segments;
		Cships[i].polygon.segments_lla.resize(n_ship_segments);
		Cships[i].polygon.segments    .resize(n_ship_segments);
	}
}

GeometricConstr::GeometricConstr(int nbrCShips, int nbrAShips)
	:
	nbrCShips( nbrCShips ),
	nbrAShips( nbrAShips ),
	Cships   ( std::vector<ship> ( nbrCShips ) ),
	Aships   ( std::vector<ship> ( nbrAShips ) ),
	n_harb_polygons( 0 ),
	n_obst_polygons( 0 ),
	n_spee_polygons( 0 ),
	n_guid_polygons( 0 ),
	lat_max(-std::numeric_limits<double>::infinity()),
	lat_min( std::numeric_limits<double>::infinity()),
	lon_max(-std::numeric_limits<double>::infinity()),
	lon_min( std::numeric_limits<double>::infinity())
{
	for(int i=0; i<nbrCShips; i++) {
		Cships[i].polygon.n_segments = n_ship_segments;
		Cships[i].polygon.segments_lla.resize(n_ship_segments);
		Cships[i].polygon.segments    .resize(n_ship_segments);
	}
	for(int i=0; i<nbrAShips; i++) {
		Aships[i].polygon.n_segments = n_ship_segments;
		Aships[i].polygon.segments_lla.resize(n_ship_segments);
		Aships[i].polygon.segments    .resize(n_ship_segments);
	}
}

void GeometricConstr::loadPolygon(std::string filename, vector<poly> *polygons, unsigned int *n_polygons){
	std::ifstream f;
	std::string str;
	f.open ( filename.c_str(), std::ios::in );
	if ( !f ) {
		throw std::runtime_error ( "Error Reading File in GeometricConstr::loadObstaclePolygon; Wrong Filename/Path?" );
	}

	std::string lon;
	std::string lat;
	getline ( f,str );
	
	while ( !f.eof() ) {
		(*n_polygons)++;
		polygons->push_back(poly());
		double center_x=0.;
		double center_y=0.;
		while( str.size()>0 ){
			unsigned int i=0;
			while ( i<str.size() && str.at ( i ) ==' ' ) {
				i++;
			}
			lon  = "";
			lat  = "";
			while ( i<str.size() && str.at ( i ) !=' ' ) {
				lon += str.at ( i );
				i++;
			}
			while ( i<str.size() && str.at ( i ) ==' ' ) {
				i++;
			}
			while ( i<str.size() && str.at ( i ) !=' ' ) {
				lat += str.at ( i );
				i++;
			}
			double lat_d = atof(lat.c_str());
			double lon_d = atof(lon.c_str());
			if(lat_max<lat_d) lat_max=lat_d;
			if(lat_min>lat_d) lat_min=lat_d;
			if(lon_max<lon_d) lon_max=lon_d;
			if(lon_min>lon_d) lon_min=lon_d;
			polygons->back().segments_lla.push_back ( Vec2d( lat_d, lon_d ) );
			polygons->back().n_segments++;
			center_x += lat_d;
			center_y += lon_d;
			getline ( f,str );
		}
		center_x /= polygons->back().n_segments;
		center_y /= polygons->back().n_segments;
		polygons->back().center_lla = Vec2d(center_x, center_y);
		getline ( f,str );
	}
	f.close();
}

void GeometricConstr::transformPolygonsToNED(vector<poly> *polygons, unsigned int *n_polygons){
	
	double z;

	for (unsigned int i=0; i<*n_polygons; i++) {
		(*polygons)[i].segments.resize((*polygons)[i].segments_lla.size());
		for (unsigned int j=0; j<(*polygons)[i].n_segments; j++) {
			LLAtoNED((*polygons)[i].segments_lla[j].x, (*polygons)[i].segments_lla[j].y, 0., lat0, lon0, 0.,
							  &((*polygons)[i].segments[j].x), &((*polygons)[i].segments[j].y), &z);
		}
		LLAtoNED((*polygons)[i].center_lla.x, (*polygons)[i].center_lla.y, 0., lat0, lon0, 0.,
						  &((*polygons)[i].center.x), &((*polygons)[i].center.y), &z);
		(*polygons)[i].area = areaOfPoly((*polygons)[i].segments);
	}
}

void GeometricConstr::initializeShip(ship *s, double length, double width) {
	s->length = length;
	s->width  = width;
	s->polygon.area = length*width;

	s->polygon.center.x=0.;
	s->polygon.center.y=0.;
	s->heading =0.;
	s->polygon.segments[0].x=0.;
	s->polygon.segments[0].y=0.;
	s->polygon.segments[1].x=0.;
	s->polygon.segments[1].y=0.;
	s->polygon.segments[2].x=0.;
	s->polygon.segments[2].y=0.;
	s->polygon.segments[3].x=0.;
	s->polygon.segments[3].y=0.;
}

void GeometricConstr::updateShip(ship *s, double center_x, double center_y, double heading) {
	s->polygon.center.x=center_x;
	s->polygon.center.y=center_y;
	s->heading =heading;

	double l_x= cos( s->heading ) * s->length/2;
	double w_x=-sin( s->heading ) * s->width/2;
	double l_y= sin( s->heading ) * s->length/2;
	double w_y= cos( s->heading ) * s->width/2;
	
	double l_x_b= cos( s->heading ) * (s->length+0*25)/2;
	double w_x_b=-sin( s->heading ) * (s->width +0*25)/2;
	double l_y_b= sin( s->heading ) * (s->length+0*25)/2;
	double w_y_b= cos( s->heading ) * (s->width +0*25)/2;

	s->polygon.segments[0].x=s->polygon.center.x +l_x+w_x;
	s->polygon.segments[0].y=s->polygon.center.y +l_y+w_y;
	s->polygon.segments[1].x=s->polygon.center.x -l_x+w_x;
	s->polygon.segments[1].y=s->polygon.center.y -l_y+w_y;
	s->polygon.segments[2].x=s->polygon.center.x -l_x-w_x;
	s->polygon.segments[2].y=s->polygon.center.y -l_y-w_y;
	s->polygon.segments[3].x=s->polygon.center.x +l_x-w_x;
	s->polygon.segments[3].y=s->polygon.center.y +l_y-w_y;
}

void GeometricConstr::loadHarbPolygon(std::string filename){
	loadPolygon(filename, &harb_polygons, &n_harb_polygons);
	harb_polygons.back().map_grid = MapGrid();
}

void GeometricConstr::loadHarbPolygon(std::string filename, std::string gridfilename){
	loadPolygon(filename, &harb_polygons, &n_harb_polygons);
	if(gridfilename == "") harb_polygons.back().map_grid = MapGrid();
	if(gridfilename != "") harb_polygons.back().map_grid = MapGrid(gridfilename, lat0, lon0, hgt0);
}

void GeometricConstr::loadObstPolygon(std::string filename){
	loadPolygon(filename, &obst_polygons, &n_obst_polygons);
	obst_polygons.back().map_grid = MapGrid();
}

void GeometricConstr::loadObstPolygon(std::string filename, std::string gridfilename){
	loadPolygon(filename, &obst_polygons, &n_obst_polygons);
	if(gridfilename == "") obst_polygons.back().map_grid = MapGrid();
	if(gridfilename != "") obst_polygons.back().map_grid = MapGrid(gridfilename, lat0, lon0, hgt0);
}

void GeometricConstr::loadSpeePolygon(std::string filename){
	if(!filename.empty()) loadPolygon(filename, &spee_polygons, &n_spee_polygons);
}

void GeometricConstr::loadGuidPolygon(std::string filename){
	if(!filename.empty())  loadPolygon(filename, &guid_polygons, &n_guid_polygons);
}

void GeometricConstr::transformHarbPolygonsToNED(){
	transformPolygonsToNED(&harb_polygons, &n_harb_polygons);
}

void GeometricConstr::transformObstPolygonsToNED(){
	transformPolygonsToNED(&obst_polygons, &n_obst_polygons);
}

void GeometricConstr::transformSpeePolygonsToNED(){
	transformPolygonsToNED(&spee_polygons, &n_spee_polygons);
}

void GeometricConstr::transformGuidPolygonsToNED(){
	transformPolygonsToNED(&guid_polygons, &n_guid_polygons);
}

void GeometricConstr::initializeCShip(int i, double length, double width) {
	initializeShip(&(Cships[i]), length, width);
}

void GeometricConstr::initializeAShip(int i, double length, double width) {
	initializeShip(&(Aships[i]), length, width);
}

void GeometricConstr::updateCShip(int i, double center_x, double center_y, double heading) {
	updateShip( &(Cships[i]), center_x, center_y, heading );
}

void GeometricConstr::updateAShip(int i, double center_x, double center_y, double heading) {
	updateShip( &(Aships[i]), center_x, center_y, heading );
}

void GeometricConstr::DistArea( poly *p1, poly *p2, double* g1, double* g2 ) {
	/*
		Calculation of the distance and intersection area between two polygons.
		g1 (>=min_dist) represents the distance between both polygons
						(negative, if one polygon is completely inside the other)
		g2 (=0) represents the intersection area between both polygons
	*/
	*g1=sqDinstancePolyPoly( p1->segments, p2->segments );
	
	// The polygons intersect eachother
	if( *g1<EPSILON ){
		std::vector<Vec2d> res;
		poly_clip( p1->segments, p2->segments, res );
		*g2=areaOfPoly( res );
		*g2=A_factor * (*g2) * (*g2);
	// The polygons do NOT intersect eachother
	}else{
		// One polygon lies completely inside the other
		if( isPointInsidePolygon(p1->segments, p2->center.x, p2->center.y ) || 
			isPointInsidePolygon(p2->segments, p1->center.x, p1->center.y )   ){
			*g1= (-1)*(*g1);
			*g2 = min2(p1->area, p2->area);
			*g2=A_factor * (*g2) * (*g2);
		// The polygons have a distance greater than zero
		}else{
			*g2=0.;
		}
	}
}

void GeometricConstr::DistArea( poly *p1, poly *p2, MapGrid *map_grid, double* g1, double* g2){
	double lat;
	double lon;
	double hgt;
	NEDtoLLA(p1->center.x, p1->center.y, 0, lat0, lon0, hgt0, &lat, &lon, &hgt);
	
	// finding the grids the ship lays in
	unsigned int i_x = (unsigned int) floor( (lat - map_grid->min_lat)/map_grid->width );
	unsigned int j_y = (unsigned int) floor( (lon - map_grid->min_lon)/map_grid->width );
	
	if(i_x>map_grid->lat_size-1){
		i_x=map_grid->lat_size-1;
	}
	if(j_y>map_grid->lon_size-1){
		j_y=map_grid->lon_size-1;
	}
	
	unsigned int min_idx_x;
	unsigned int max_idx_x;
	unsigned int min_idx_y;
	unsigned int max_idx_y;
	
	if( fmod((p1->center.x - map_grid->min_lat), map_grid->width)<map_grid->width/2.){
		if(i_x>0){
			min_idx_x=i_x-1;
			max_idx_x=i_x;
		}else{
			min_idx_x=i_x;
			max_idx_x=i_x;
		}
	}else{
		if(i_x<map_grid->lat_size-1){
			min_idx_x=i_x;
			max_idx_x=i_x+1;
		}else{
			min_idx_x=i_x;
			max_idx_x=i_x;
		}
	}
	if( fmod((p1->center.y - map_grid->min_lon), map_grid->width)<map_grid->width/2.){
		if(j_y>0){
			min_idx_y=j_y-1;
			max_idx_y=j_y;
		}else{
			min_idx_y=j_y;
			max_idx_y=j_y;
		}
	}else{
		if(j_y<map_grid->lon_size-1){
			min_idx_y=j_y;
			max_idx_y=j_y+1;
		}else{
			min_idx_y=j_y;
			max_idx_y=j_y;
		}
	}
	
	// calculation of distance
	double g1_aux=0.;
	*g1=INFINITY;
	for (unsigned int i=min_idx_x; i<=max_idx_x; i++) {
		for (unsigned int j=min_idx_y; j<=max_idx_y; j++) {
			for (unsigned int m=0; m<map_grid->segs_in_dist[i][j].size(); m++) {
				g1_aux = sqDinstancePolySegment( p1->segments, p2->segments , map_grid->segs_in_dist[i][j][m] );
				if(g1_aux < *g1){
					*g1 = g1_aux;
				}
			}
		}
	}
	
	// The polygons intersect eachother
	if( *g1<EPSILON ){
		*g2=0.;
		std::vector<Vec2d> res;
		for (unsigned int i=min_idx_x; i<=max_idx_x; i++) {
			for (unsigned int j=min_idx_y; j<=max_idx_y; j++) {
				if( map_grid->harbor_area[i][j].size()>0 ){
					res.clear();
					poly_clip( p1->segments, map_grid->harbor_area[i][j], res );
					*g2 += areaOfPoly( res );
				}
			}
		}
		*g2=A_factor * (*g2) * (*g2);
// 			*g2=A_factor*(*g2)*tanh(*g2);
		
	// The polygons do NOT intersect eachother
	}else{
// 			*g1 = sqrt(*g1);
// 			*g1 = *g1*tanh(*g1);
		bool is_inside=false;
		for (unsigned int i=min_idx_x; !is_inside && i<=max_idx_x; i++) {
			for (unsigned int j=min_idx_y; !is_inside && j<=max_idx_y; j++) {
				if( map_grid->harbor_area[i][j].size()>0 ){
					is_inside = isPointInsidePolygon( map_grid->harbor_area[i][j], p1->center.x, p1->center.y );
				}
			}
		}
		// The ship lies completely inside the harbor:
		if( is_inside ){
			*g1= (-1)*(*g1);
			*g2 = min2(p1->area, p2->area);
			*g2 = A_factor*(*g2)*(*g2);
// 				*g2 = A_factor*(*g2)*tanh(*g2);
		
		// The ship lies completely outside the harbor
		}else{
			*g2=0.;
		}
	}
}

void GeometricConstr::ConstrCShipCShip( unsigned int k, unsigned int l, double* g1, double* g2, bool withGrid ) {
	/*
		Calculation of the collision avoidance of ship k and ship l.
		g1 (>=min_dist) represents the distance between both ships
						(negative, if one ship is completely inside the other)
		g2 (=0) represents the intersection area between both ships
	*/
	DistArea( &(Cships[k].polygon), &(Cships[l].polygon), g1, g2);
}

void GeometricConstr::ConstrCShipAShip( unsigned int k, unsigned int l, double* g1, double* g2, bool withGrid ) {
	/*
		Calculation of the collision avoidance of ship k and ship l.
		g1 (>=min_dist) represents the distance between both ships
						(negative, if one ship is completely inside the other)
		g2 (=0) represents the intersection area between both ships
	*/
	DistArea( &(Cships[k].polygon), &(Aships[l].polygon), g1, g2);
}

void GeometricConstr::ConstrCShipObsta( unsigned int k, unsigned int l, double* g1, double* g2, bool withGrid ) {
	/*
		Calculation of the collision avoidance of ship k and obstacle l.
		g1 (>=min_dist) represents the distance between both objects
						(negative, if one is completely inside the other)
		g2 (=0) represents the intersection area between both objects
	*/
	if(withGrid) {
		DistArea( &(Cships[k].polygon), &(obst_polygons[l]), &(obst_polygons[l].map_grid), g1, g2);
	}
	// calculation with complete harbor-polygon
	if(!withGrid) {
		DistArea( &(Cships[k].polygon), &(obst_polygons[l]), g1, g2);
	}
}

void GeometricConstr::ConstrCShipHarbr( unsigned int k, unsigned int l, double* g1, double* g2, bool withGrid ) {
	/*
		Calculation of the collision avoidance of ship k and the harbor.
		g1 (>=min_dist) represents the distance between ship and harbor
						(negative, if the ship is completely outside the harbor)
		g2 (>=A) represents the intersection area between ship and harbor
	*/
    //calculation with harbor-polygon segmented in map-grid
    if(withGrid) {
		DistArea( &(Cships[k].polygon), &(harb_polygons[l]), &(harb_polygons[l].map_grid), g1, g2);
	}
	// calculation with complete harbor-polygon
	if(!withGrid) {
		DistArea( &(Cships[k].polygon), &(harb_polygons[l]), g1, g2);
	}
}



void GeometricConstr::sqDistShipGuide( unsigned int k, double* g1 ) {
	/*
		Calculation of the square of the distance between ship k and the guiding polygon.
	*/

	Vec2d pos(Cships[k].polygon.center.x, Cships[k].polygon.center.y);
	*g1 = sqDistancePolyPoint( guid_polygons[0].segments, pos );//has ato be fixed
}



double GeometricConstr::SpeedBound( unsigned int k,vector<double> speedbounds ) {
	double g1;
	double g2;
	double sqship_area = getAreaOfShip(k)*getAreaOfShip(k);
	double bound = 0;
	
	for(int l=0; l<n_spee_polygons; l++) {
		
		DistArea( &(Cships[k].polygon), &(spee_polygons[l]), &g1, &g2);
		bound += g2/sqship_area/A_factor*(100 - speedbounds[l]);

	}
	return 100-bound;
}


bool GeometricConstr::isPointInsidePolygon( std::vector<Vec2d>& segs, double pt_x, double pt_y ) {
	
// Copyright (c) 1970-2003, Wm. Randolph Franklin
	
	unsigned int i;
	unsigned int j;
	bool within = false;

	for ( i = 0, j = segs.size()-1; i < segs.size(); j = i++ ) {
		if ( ( ( segs[i].y>pt_y ) != ( segs[j].y>pt_y ) ) &&
				( pt_x < ( segs[j].x-segs[i].x ) * ( pt_y-segs[i].y ) / ( segs[j].y-segs[i].y ) + segs[i].x ) ) {
			within = !within;
		}
	}

	return within;
}

double GeometricConstr::sqDinstancePolyPoly( std::vector<Vec2d>& segs1, std::vector<Vec2d>& segs2 ) {

	double dist=INFINITY;
	double new_dist=0.;
	
	unsigned int n_segs1 = segs1.size();
	unsigned int n_segs2 = segs2.size();

	// Calculating the squared distance between both polygons, which corresponds
	// to the minimal distance between all pairs of segments
	for ( unsigned int i=0; i<n_segs1-1; i++ ) {
		for ( unsigned int j=0; j<n_segs2-1; j++ ) {
			new_dist=sqDinstanceSegmentSegment( segs1[i], segs1[i+1],
												segs2[j], segs2[j+1] );

			if ( new_dist<dist ) {
				dist=new_dist;
			}
		}
	}

	for ( unsigned int j=0; j<n_segs2-1; j++ ) {
		new_dist=sqDinstanceSegmentSegment( segs1[n_segs1-1],
											segs1[0],
											segs2[j], segs2[j+1] );

		if ( new_dist<dist ) {
			dist=new_dist;
		}
	}

	new_dist=sqDinstanceSegmentSegment( segs1[n_segs1-1],
										segs1[0],
										segs2[n_segs2-1],
										segs2[0] );

	if ( new_dist<dist ) {
		dist=new_dist;
	}

	return dist;
}

double GeometricConstr::sqDinstancePolySegment( std::vector<Vec2d>& segs1, std::vector<Vec2d>& segs2, unsigned int segs2_idx ) {

	double dist=INFINITY;
	double new_dist=0.;
	
	unsigned int n_segs1 = segs1.size();

	// Calculating the squared distance between polygon and segment, which corresponds
	// to the minimal distance between all polygon segments and the single segment
	for ( unsigned int i=0; i<n_segs1-1; i++ ) {
		new_dist=sqDinstanceSegmentSegment( segs1[i], segs1[i+1],
											segs2[segs2_idx], segs2[segs2_idx+1] );

		if ( new_dist<dist ) {
			dist=new_dist;
		}
	}
	new_dist=sqDinstanceSegmentSegment( segs1[n_segs1-1],
										segs1[0],
										segs2[segs2_idx], segs2[segs2_idx+1] );

	if ( new_dist<dist ) {
		dist=new_dist;
	}

	return dist;
}

double GeometricConstr::sqDinstanceSegmentSegment( Vec2d &p1, Vec2d &q1, Vec2d &p2, Vec2d &q2 ) {
	Vec2d d1 = q1 - p1; // Direction vector of segment S1
	Vec2d d2 = q2 - p2; // Direction vector of segment S2
	Vec2d r = p1 - p2;
	double a = d1.dot( d1 );  // Squared length of segment S1, always nonnegative
	double e = d2.dot( d2 );  // Squared length of segment S2, always nonnegative
	double f = d2.dot( r );
	double s=0.;
	double t=0.;

	// Check if either or both segments degenerate into points
	if ( a <= EPSILON && e <= EPSILON ) {
		// Both segments degenerate into points
		return ( p1-p2 ).dot( p1-p2 );
	}

	if ( a <= EPSILON ) {
		// First segment degenerates into a point
		s = 0.0;
		t = f / e; // s = 0 => t = (b*s + f) / e = f / e
		t = Clamp( t, 0.0, 1.0 );
	} else {
		double c = d1.dot( r );

		if ( e <= EPSILON ) {
			// Second segment degenerates into a point
			t = 0.0;
			s = Clamp( -c / a, 0.0, 1.0 );  // t = 0 => s = (b*t - c) / a = -c / a
		} else {
			// The general nondegenerate case starts here
			double b = d1.dot( d2 );
			double denom = a*e-b*b; // Always nonnegative

			// If segments not parallel, compute closest point on L1 to L2 and
			// clamp to segment S1. Else pick arbitrary s (here 0)
			if ( denom != 0.0 ) {
				s = Clamp( ( b*f - c*e ) / denom, 0.0, 1.0 );
			} else s = 0.0;

			// Compute point on L2 closest to S1(s) using
			// t = Dot((P1 + D1*s) - P2,D2) / Dot(D2,D2) = (b*s + f) / e
			double tnom = b*s + f;

			if ( tnom < 0.0 ) {
				t = 0.0;
				s = Clamp( -c / a, 0.0, 1.0 );
			} else if ( tnom > e ) {
				t = 1.0;
				s = Clamp( ( b - c ) / a, 0.0, 1.0 );
			} else {
				t = tnom / e;
			}
		}
	}

	return ( ( p1 + d1 * s )- ( p2 + d2 * t ) ).dot( ( p1 + d1 * s )- ( p2 + d2 * t ) );
}

double GeometricConstr::sqDistancePolyPoint( std::vector<Vec2d>& segs1, Vec2d &p2 ) {

	double dist=INFINITY;
	double new_dist=0.;
	
	unsigned int n_segs1 = segs1.size();

	// Calculating the squared distance between polygon and point, which corresponds
	// to the minimal distance between all polygon segments and the single point
	for ( unsigned int i=0; i<n_segs1-1; i++ ) {
		new_dist=sqDinstanceSegmentPoint( segs1[i], segs1[i+1], p2 );

		if ( new_dist<dist ) {
			dist=new_dist;
		}
	}
	new_dist=sqDinstanceSegmentPoint( segs1[n_segs1-1], segs1[0], p2 );

	if ( new_dist<dist ) {
		dist=new_dist;
	}

	return dist;
}

double GeometricConstr::sqDinstanceSegmentPoint( Vec2d &p1, Vec2d &q1, Vec2d &p2 ) {
	Vec2d d1 = q1 - p1; // Direction vector of segment S1
	Vec2d r = p1 - p2;
	double a = d1.dot( d1 );  // Squared length of segment S1, always nonnegative
	double s=0.;
	
	// Check if segment degenerates into a point
	if ( a <= EPSILON ) {
		// Segment degenerates into a point
		return ( p1-p2 ).dot( p1-p2 );
	}

	double c = d1.dot( r );

	s = Clamp( -c / a, 0.0, 1.0 );  // t = 0 => s = (b*t - c) / a = -c / a

	return ( ( p1 + d1 * s ) - p2  ).dot( ( p1 + d1 * s ) - p2 );
}

double GeometricConstr::Clamp( double n, double min, double max ) {
	if ( n < min ) return min;
	if ( n > max ) return max;
	return n;
}

double GeometricConstr::areaOfPoly( std::vector<Vec2d>& segs ) {
	double A=0.;

	if ( segs.size()>2 ) {
		for ( unsigned int i=0; i<segs.size()-1; i++ ) {
			A=A+( ( segs[i].x * segs[i+1].y ) - ( segs[i].y * segs[i+1].x ) ); 
		}
		A=A+( ( segs[segs.size()-1].x * segs[0].y ) - ( segs[segs.size()-1].y * segs[0].x ) ); 
		A=fabs( A )/2.;
	}
	return A;
}

int GeometricConstr::line_sect( Vec2d& x0, Vec2d& x1, Vec2d& y0, Vec2d& y1, Vec2d& res ) {
	Vec2d dx=x1-x0;
	Vec2d dy=y1-y0;
	Vec2d d=x0-y0;
	/* x0 + a dx = y0 + b dy ->
	   x0 X dx = y0 X dx + b dy X dx ->
	   b = (x0 - y0) X dx / (dy X dx) */
	double dyx = dy.cross( dx );

	if ( !dyx ) return 0;

	dyx = d.cross( dx ) / dyx;

	if ( dyx <= 0 || dyx >= 1 ) return 0;

	res.x = y0.x + dyx * dy.x;
	res.y = y0.y + dyx * dy.y;
	return 1;
}

/* tells if vec c lies on the left side of directed edge a->b
 * 1 if left, -1 if right, 0 if colinear
 */
int GeometricConstr::left_of( Vec2d& a, Vec2d& b, Vec2d& c ) {
	double out = ( b-a ).cross( c-b );
	return out < 0 ? -1 : out > 0;
}



/* this works only if all of the following are true:
 *   1. poly has no colinear edges;
 *   2. poly has no duplicate vertices;
 *   3. poly has at least three vertices;
 *   4. poly is convex (implying 3).
*/
int GeometricConstr::poly_winding( std::vector<Vec2d>& p ) {
	return left_of( p[0], p[1], p[2] );
}

void GeometricConstr::poly_edge_clip( std::vector<Vec2d>& sub, Vec2d& x0, Vec2d& x1, int left, std::vector<Vec2d>& res ) {
	int i, side0, side1;
	Vec2d tmp;
	Vec2d v0 = sub.back();
	res.resize( 0 );

	side0 = left_of( x0, x1, v0 );

	if ( side0 != -left ) res.push_back( v0 );

	for ( i = 0; i < sub.size(); i++ ) {
		side1 = left_of( x0, x1, sub[i] );

		if ( side0 + side1 == 0 && side0 )

			/* last point and current straddle the edge */
			if ( line_sect( x0, x1, v0, sub[i], tmp ) )
				res.push_back( tmp );

		if ( i == sub.size()-1 ) break;

		if ( side1 != -left ) res.push_back( sub[i] );

		v0 = sub[i];
		side0 = side1;
	}
}

void GeometricConstr::poly_clip( std::vector<Vec2d>& sub, std::vector<Vec2d>& clip, std::vector<Vec2d>& p2) {
	std::vector<Vec2d> p1;
	p2.resize(0);
	std::vector<Vec2d> tmp;

	int dir = poly_winding( clip );
	poly_edge_clip( sub, clip.back(), clip[0], dir, p2 );

	for ( int i = 0; i < clip.size()-1; i++ ) {
		tmp = p2;
		p2 = p1;
		p1 = tmp;

		if ( p1.size() == 0 ) {
			p2.resize( 0 );
			break;
		}

		poly_edge_clip( p1, clip[i], clip[i+1], dir, p2 );
	}
}

std::vector<Vec2d>* GeometricConstr::getHarborSegements(unsigned int l){
	return &(harb_polygons[l].segments);
}

std::vector<Vec2d>* GeometricConstr::getObstaclesSegements(unsigned int l){
	return &(obst_polygons[l].segments);
}

unsigned int GeometricConstr::getNbrObstacles(){
	return n_obst_polygons;
}
		
double GeometricConstr::getAFactor(){
	return A_factor;
}

double GeometricConstr::getAreaOfShip(unsigned int i){
	return Cships[i].polygon.area;
}

void GeometricConstr::initializeOccMap(Eigen::Vector2d cell_dim_ned, double buffer){
	if(n_harb_polygons<=0){
		assert("Cannot initialize occupancy map, no harbor-polygon loaded");
	}

	occ_map_cell_dim_ned=cell_dim_ned;

	double x_l_ned = 0.;
	double x_u_ned = 0.;
	double y_l_ned = 0.;
	double y_u_ned = 0.;
	double dummy = 0.;
	
	LLAtoNED(lat_min, lon_min, 0., lat0, lon0, hgt0,
			&x_l_ned, &y_l_ned, &dummy);
	LLAtoNED(lat_max, lon_max, 0., lat0, lon0, hgt0,
			&x_u_ned, &y_u_ned, &dummy);
	dsl::Vector2i nr_cells( floor( (x_u_ned-x_l_ned)/cell_dim_ned[0] ),
							floor( (y_u_ned-y_l_ned)/cell_dim_ned[1] ) );

	// cout<<setprecision(10)
	// 	<<"x="<<x_l_ned<<" - "<<x_u_ned<<endl
	// 	<<"y="<<y_l_ned<<" - "<<y_u_ned<<endl<<endl;

	if(create_new_occ_map){
		occ_map = new dsl::TypedMap ( Eigen::Vector2d(x_l_ned, y_l_ned),
									Eigen::Vector2d(x_u_ned, y_u_ned), nr_cells );

		for(int id=0; id<occ_map->nc; id++){
			std::vector<Eigen::Vector2d> pts_eigen = occ_map->Pos ( id );

			// conversion from Eigen::Vector2d-ned to Vec2d-ned
			std::vector<Vec2d> pts_ned(pts_eigen.size());
			double center_x=0.;
			double center_y=0.;
			double dummy=0.;
			for(int j=0; j<pts_eigen.size(); j++){
				pts_ned[j].x=pts_eigen[j][0];
				pts_ned[j].y=pts_eigen[j][1];
				center_x += pts_ned[j].x;
				center_y += pts_ned[j].y;
			}
			center_x /= pts_ned.size();
			center_y /= pts_ned.size();
			bool occupied = false;
			double dist = sqDinstancePolyPoly(pts_ned, harb_polygons[0].segments);
			bool is_inside = isPointInsidePolygon(harb_polygons[0].segments, center_x, center_y);
			if( dist<buffer*buffer+EPSILON || !is_inside) occupied=true;
			for(int j=0; j<n_obst_polygons; j++){
				dist = sqDinstancePolyPoly(pts_ned, obst_polygons[j].segments);
				is_inside = isPointInsidePolygon(obst_polygons[j].segments, center_x, center_y);
				if( dist<buffer*buffer+EPSILON || is_inside ) occupied=true;
			}
			occ_map->Set(id, occupied);
		}

		occ_map->empty = true;
		
		// save occupancy map as image
		string filename = "occ_map.ppm";
		save(*occ_map, filename.c_str() );
		
	}else{

		double x_l_ned = 0.;
		double x_u_ned = 0.;
		double y_l_ned = 0.;
		double y_u_ned = 0.;
		double dummy = 0.;
		LLAtoNED(lat_min, lon_min, 0., lat0, lon0, hgt0,
				&x_l_ned, &y_l_ned, &dummy);
		LLAtoNED(lat_max, lon_max, 0., lat0, lon0, hgt0,
				&x_u_ned, &y_u_ned, &dummy);

		Eigen::Vector2d cell_dim_ned;
		cell_dim_ned[0]=(x_u_ned-x_l_ned)/nr_cells[0];
		cell_dim_ned[1]=(y_u_ned-y_l_ned)/nr_cells[1];

		dsl::TypedMap map_tmp = dsl::load ( "occ_map.ppm", cell_dim_ned );
		occ_map = new dsl::TypedMap ( map_tmp );

		occ_map->empty = true;
				
		occ_map->xlb=Eigen::Vector2d(x_l_ned, y_l_ned);
		occ_map->xub=Eigen::Vector2d(x_u_ned, y_u_ned);
	}

	// cout<<endl<<lat_min<<"   "<<lat_max<<endl;
	// cout      <<lon_min<<"   "<<lon_max<<endl;
}









