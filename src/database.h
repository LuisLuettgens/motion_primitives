#ifndef DATABASE_H
#define DATABASE_H

#include <deque>	// like a double-sided vector. You can perform pushfront(), pushback(), popfront(), popback()
#include <vector>
#include "geometric_constr.h"
#include "Models/BaseModel.h"

enum SHIP_ID {
	MODEL1V1,
	MODEL1V2,
	MODEL1V3,
	MODEL1V4,
	MODEL1V5,
	MODEL1V5MECKPOMM,
	MODEL2V1,
	MODELAIS
};

using namespace std;

struct ship {
	string id;
	BaseModel* model;
};

struct polygon {
	
};

struct COOP {
	deque<ship> ships_considered;
	deque<string> ships_to_include;
	deque<string> ships_to_remove;
};

struct AIS : COOP {				// AIS inherits from COOP and adds additional members
	double time_stamp;
	double lat;
	double lon;
	double sog;
	double cog;
	vector<Vec2d> waypoints;
};

struct database {
	deque<ship> ship_dictionary;		// here all ships shall be stored
	COOP coop;				// these are the cooperative ships
	AIS ais;				// these are the non-cooperative ships

	vector<Vec2d> harbour_polygon_lla;
};


struct ship2 {
	SHIP_ID id;
	BaseModel* model;
};

struct database2 {
	deque<ship> ship_dictionary;
	
	deque<ship2> ships_considered;
	deque<SHIP_ID> ships_to_include;
	deque<SHIP_ID> ships_to_remove;

	deque<ship2> AIS_considered;
	deque<SHIP_ID> AIS_to_include;
	deque<SHIP_ID> AIS_to_remove;

	vector<Vec2d> harbour_polygon_lla;
};

#endif