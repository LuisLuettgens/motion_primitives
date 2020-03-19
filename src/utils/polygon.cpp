#include "polygon.h"
// #include "../global.h"


Poly::Poly():nbrVertices(0){}
	
void Poly::LLAtoNED_(double lat0, double lon0, double hgt0){
	double z;
	for (unsigned int i=0; i<nbrVertices; i++) {
		LLAtoNED(vertices_lla[i].x, vertices_lla[i].y, 0., lat0, lon0, hgt0, &(vertices_ned[i].x), &(vertices_ned[i].y), &z);
	}
	LLAtoNED(center_lla.x, center_lla.y, 0., lat0, lon0, 0., &(center_ned.x), &(center_ned.y), &z);
}
	
void Poly::NEDtoLLA_(double lat0, double lon0, double hgt0){
	double z;
	for (unsigned int i=0; i<nbrVertices; i++) {
		NEDtoLLA(vertices_ned[i].x, vertices_ned[i].y, 0., lat0, lon0, hgt0, &(vertices_lla[i].x), &(vertices_lla[i].y), &z);
	}
	NEDtoLLA(center_ned.x, center_ned.y, 0., lat0, lon0, 0., &(center_lla.x), &(center_lla.y), &z);
}
	
double D(Poly p){
	return sqrt(sqD(p));
}
	
double Poly::sqD(Poly p){
	double dist=INFINITY;
	double new_dist=0.;

	// Calculating the squared distance between both polygons, which corresponds
	// to the minimal distance between all pairs of segments
	for ( unsigned int i=0; i<nbrVertices-1; i++ ) {
		for ( unsigned int j=0; j<p.nbrVertices-1; j++ ) {
			new_dist=sqDinstanceSegmentSegment(   vertices_ned[i],   vertices_ned[i+1],
												p.vertices_ned[j], p.vertices_ned[j+1] );

			if ( new_dist<dist ) {
				dist=new_dist;
			}
		}
	}

	for ( unsigned int j=0; j<p.nbrVertices-1; j++ ) {
		new_dist=sqDinstanceSegmentSegment(   vertices_ned[nbrVertices-1],   vertices_ned[0],
											p.vertices_ned[j]            , p.vertices_ned[j+1] );

		if ( new_dist<dist ) {
			dist=new_dist;
		}
	}

	new_dist=sqDinstanceSegmentSegment(   vertices_ned[  nbrVertices-1],   vertices_ned[0],
										p.vertices_ned[p.nbrVertices-1], p.vertices_ned[0] );

	if ( new_dist<dist ) {
		dist=new_dist;
	}

	return dist;
}
	
double Poly::Area_( ) {
	double area=0.;
	if ( nbrVertices>2 ) {
		for ( unsigned int i=0; i<nbrVertices-1; i++ ) {
			area=area+( ( vertices_ned[i].x * vertices_ned[i+1].y ) - ( vertices_ned[i].y * vertices_ned[i+1].x ) ); 
		}
		area=area+( ( vertices_ned[nbrVertices-1].x * vertices_ned[0].y ) - ( vertices_ned[nbrVertices-1].y * vertices_ned[0].x ) );
		area=fabs( area )/2.;
	}
	return area;
}
	
double Poly::interA(Poly poly){
	double sqdist=sqD(poly);

	// The polygons intersect eachother
	if( sqdist<EPSILON ){
		std::vector<Vec2d> res;
		poly_clip( poly.vertices_ned, vertices_ned, res );
		return Area( res );
	// The polygons do NOT intersect eachother
	}else{
		// One object lies completely inside the other
		if( isPointInside( poly.center_ned ) || poly.isPointInside( center_ned ) ){
			std::vector<Vec2d> res;
			poly_clip( poly.vertices_ned, vertices_ned, res );
			return Area( res );
		// The objects have a positive distance greater than zero
		}else{
			return 0.;
		}
	}
}
	
bool Poly::isPointInside( Vec2d pt ) {
// Copyright (c) 1970-2003, Wm. Randolph Franklin
	
	unsigned int i;
	unsigned int j;
	bool within = false;

	for ( i = 0, j = nbrVertices-1; i < nbrVertices; j = i++ ) {
		if ( ( ( vertices_ned[i].y>pt.y ) != ( vertices_ned[j].y>pt.y ) ) &&
				( pt.x < ( vertices_ned[j].x-vertices_ned[i].x ) * ( pt.y-vertices_ned[i].y ) / ( vertices_ned[j].y-vertices_ned[i].y ) + vertices_ned[i].x ) ) {
			within = !within;
		}
	}
	return within;
}

