#ifndef functions_h
#define functions_h

#include "TransWORHP.h"
#include "TWproblem.h"
#include "TWfolder.h"

#include "../Models/Model1v5MeckPomm.h"
#include "../Models/Buggy_v1.h"


using namespace std;

void convert(string datatype, int *v);

struct item
{
	string Description;
	string Name;
	int value;
};

struct Description
{
	int    ID;
	string Datatype;
	double Min;
	double Max;
};

struct Range
{
	int Size;
	long int IntMin;
	long int IntMax;
};

struct AISShip
{
	int ID;
	double x0;
	double y0;
	double psi;
	double cog;
	double sog;
	double r;
	double accuracy;
	double integrity;
};

struct AISShipDim
{
	double L;
	double BOA;
};

double min2(double a, double b);

class Vec2d {
	public:
		Vec2d()
		: x(0.), y(0.)
		{}
		Vec2d(double x_, double y_)
		: x(x_), y(y_)
		{}
		double x;
		double y;
		double dot(Vec2d in){
			return (x*in.x+y*in.y);
		}
		double cross(Vec2d in){
			return (x*in.y-y*in.x);
		}
		Vec2d operator - (Vec2d right){
			return Vec2d(x-right.x,y-right.y);
		}
		Vec2d operator + (Vec2d right){
			return Vec2d(x+right.x,y+right.y);
		}
		Vec2d operator * (double fac){
			return Vec2d(x*fac,y*fac);
		}
		friend std::ostream& operator<<(std::ostream& os, const Vec2d& vec){
			os << vec.x << " " << vec.y;
			return os;
		}
};



void LLAtoNED(double lat, double lon, double h, double lat0, double lon0, double h0,
			  double *x, double *y, double *z);

void NEDtoLLA(double x, double y, double z, double lat0, double lon0, double h0, 
			  double *lat, double *lon, double *h);


typedef enum {
	DP_0XA0     = 160,
	DP_0XA1     = 161,
	DP_0XA2     = 162,
	DP_0XA3     = 163,
	DP_0XA4     = 164,
	DP_0XAC     = 172,
	DP_0XAD     = 173,
	DP_0XAE     = 174,
	DP_0XAF     = 175,
} DP_ID;

typedef enum {
	NORMAL      = 0,
	SINGLEPOINT = 1,
	EXTRAPOLATE = 2
} CubicSplineBool;

typedef enum {
	COMMON_AND_FIXED       = 0,
	OPEN_WITH_FIXED_RATIOS = 1,
	OPEN                   = 2
} EndTimeBool;

typedef enum {
	FROM_MATLAB = 0,
	FROM_INTEGRATOR = 1
} InitGuessBool;

typedef enum {
	SHIP  = 0,
	BUGGY = 1
} TypeOfVehicle;

// sin in degrees
double sind(double x);

//cos in degrees
double cosd(double x);



unsigned nChoosek ( unsigned n, unsigned k );
// 
void cubicSpline(int i, int k, double dt, double dtfull, bool extrapolate, double *x0, double *y0, double *psi);
void cubicSpline(int i, int k, double dt, double dtfull, bool extrapolate, double *x0, double *y0, double *psi, double *delta, double *n);
void cubicSpline_structure(tw::DiffStructure &s, int i, int k, int idx, CubicSplineBool disc);

double t_max();
void SetFloatTime(tw::Viewer *viewer, double t_scale);

BaseModel* get_ship_model(string ID);
BaseModel* get_buggy_model(string ID);

double scaleplot (int &len, int &ndgl, int &nsteuer, double *t,double *x, double *u, int &j, int &index);

#endif




// const double EPSILON = 1.0E-5;
// double Clamp( double n, double min, double max );
// int line_sect(Vec2d& x0, Vec2d& x1, Vec2d& y0, Vec2d& y1, Vec2d& res);
// int left_of(Vec2d& a, Vec2d& b, Vec2d& c);
// double sqDinstanceSegmentSegment( Vec2d &p1, Vec2d &q1, Vec2d &p2, Vec2d &q2 );
// int    poly_winding(std::vector<Vec2d>& p);
// void   poly_edge_clip(std::vector<Vec2d>& sub, Vec2d& x0, Vec2d& x1, int left, std::vector<Vec2d>& res);
// void   poly_clip(std::vector<Vec2d>& sub, std::vector<Vec2d>& clip, std::vector<Vec2d>& p2);
// double areaOfPoly(std::vector<Vec2d>& segs);
// 
// 
// class Poly{
// public:
// 	Poly():nbrVertices(0){}
// 	unsigned int       nbrVertices;
// 	Vec2d              center_lla;
// 	Vec2d              center_ned;
// 	std::vector<Vec2d> vertices_lla;
// 	std::vector<Vec2d> vertices_ned;
// 	
// 	void LLAtoNED_(double lat0, double lon0, double hgt0){
// 		double z;
// 		for (unsigned int i=0; i<nbrVertices; i++) {
// 			LLAtoNED(vertices_lla[i].x, vertices_lla[i].y, 0., lat0, lon0, hgt0, &(vertices_ned[i].x), &(vertices_ned[i].y), &z);
// 		}
// 		LLAtoNED(center_lla.x, center_lla.y, 0., lat0, lon0, 0., &(center_ned.x), &(center_ned.y), &z);
// 	}
// 	
// 	void NEDtoLLA_(double lat0, double lon0, double hgt0){
// 		double z;
// 		for (unsigned int i=0; i<nbrVertices; i++) {
// 			NEDtoLLA(vertices_ned[i].x, vertices_ned[i].y, 0., lat0, lon0, hgt0, &(vertices_lla[i].x), &(vertices_lla[i].y), &z);
// 		}
// 		NEDtoLLA(center_ned.x, center_ned.y, 0., lat0, lon0, 0., &(center_lla.x), &(center_lla.y), &z);
// 	}
// 	
// 	double D(Poly p){
// 		return sqrt(sqD(p));
// 	}
// 	
// 	double sqD(Poly p){
// 		double dist=INFINITY;
// 		double new_dist=0.;
// 
// 		// Calculating the squared distance between both polygons, which corresponds
// 		// to the minimal distance between all pairs of segments
// 		for ( unsigned int i=0; i<nbrVertices-1; i++ ) {
// 			for ( unsigned int j=0; j<p.nbrVertices-1; j++ ) {
// 				new_dist=sqDinstanceSegmentSegment(   vertices_ned[i],   vertices_ned[i+1],
// 													p.vertices_ned[j], p.vertices_ned[j+1] );
// 
// 				if ( new_dist<dist ) {
// 					dist=new_dist;
// 				}
// 			}
// 		}
// 
// 		for ( unsigned int j=0; j<p.nbrVertices-1; j++ ) {
// 			new_dist=sqDinstanceSegmentSegment(   vertices_ned[nbrVertices-1],   vertices_ned[0],
// 												p.vertices_ned[j]            , p.vertices_ned[j+1] );
// 
// 			if ( new_dist<dist ) {
// 				dist=new_dist;
// 			}
// 		}
// 
// 		new_dist=sqDinstanceSegmentSegment(   vertices_ned[  nbrVertices-1],   vertices_ned[0],
// 											p.vertices_ned[p.nbrVertices-1], p.vertices_ned[0] );
// 
// 		if ( new_dist<dist ) {
// 			dist=new_dist;
// 		}
// 
// 		return dist;
// 	}
// 	
// 	double Area( ) {
// 		double area=0.;
// 		if ( nbrVertices>2 ) {
// 			for ( unsigned int i=0; i<nbrVertices-1; i++ ) {
// 				area=area+( ( vertices_ned[i].x * vertices_ned[i+1].y ) - ( vertices_ned[i].y * vertices_ned[i+1].x ) ); 
// 			}
// 			area=area+( ( vertices_ned[nbrVertices-1].x * vertices_ned[0].y ) - ( vertices_ned[nbrVertices-1].y * vertices_ned[0].x ) );
// 			area=fabs( area )/2.;
// 		}
// 		return area;
// 	}
// 	
// 	double interA(Poly poly){
// // 		double sqdist=sqD(poly);
// 		
// 		std::vector<Vec2d> res;
// 		poly_clip( poly.vertices_ned, vertices_ned, res );
// 		return areaOfPoly( res );
// 
// // 		// The polygons intersect eachother
// // 		if( sqdist<EPSILON ){
// // 			
// // 			poly_clip( poly.vertices_ned, vertices_ned, res );
// // 			return areaOfPoly( res );
// // 		// The polygons do NOT intersect eachother
// // 		}else{
// // 			// One object lies completely inside the other
// // 			if( isPointInside( poly.center_ned ) || poly.isPointInside( center_ned ) ){
// // 				std::vector<Vec2d> res;
// // 				poly_clip( poly.vertices_ned, vertices_ned, res );
// // 				return areaOfPoly( res );
// // 			// The objects have a positive distance greater than zero
// // 			}else{
// // 				return 0.;
// // 			}
// // 		}
// 	}
// 	
// 	bool isPointInside( Vec2d pt ) {
// 	// Copyright (c) 1970-2003, Wm. Randolph Franklin
// 		
// 		unsigned int i;
// 		unsigned int j;
// 		bool within = false;
// 
// 		for ( i = 0, j = nbrVertices-1; i < nbrVertices; j = i++ ) {
// 			if ( ( ( vertices_ned[i].y>pt.y ) != ( vertices_ned[j].y>pt.y ) ) &&
// 					( pt.x < ( vertices_ned[j].x-vertices_ned[i].x ) * ( pt.y-vertices_ned[i].y ) / ( vertices_ned[j].y-vertices_ned[i].y ) + vertices_ned[i].x ) ) {
// 				within = !within;
// 			}
// 		}
// 		return within;
// 	}
// };

// double Area( Poly poly );



