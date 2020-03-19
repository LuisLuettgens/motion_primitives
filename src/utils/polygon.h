#ifndef polygon_h
#define polygon_h

// #include "TransWORHP.h"
// #include "TWproblem.h"
// #include "TWfolder.h"
#include "functions.h"

#include "../Models/Model1v5MeckPomm.h"


using namespace std;

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


class Poly{
public:
	Poly():nbrVertices(0){}
	unsigned int       nbrVertices;
	Vec2d              center_lla;
	Vec2d              center_ned;
	std::vector<Vec2d> vertices_lla;
	std::vector<Vec2d> vertices_ned;
	
	void LLAtoNED_(double lat0, double lon0, double hgt0);
	void NEDtoLLA_(double lat0, double lon0, double hgt0);
	double D(Poly p);
	double sqD(Poly p);
	double Area_( );
	double interA(Poly poly);
	bool isPointInside( Vec2d pt );
};


#endif
