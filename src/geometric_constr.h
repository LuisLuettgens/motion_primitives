#ifndef geometric_constr_h
#define geometric_constr_h

#include "mapgrid_lla.h"
#include "utils/functions.h"
#include "data_types.h"
#include "map.h"
#include <Eigen/Geometry>
#include "utils.h"


class GeometricConstr {
	public:
		struct poly{
			poly():n_segments(0){}
			Vec2d center_lla;
			Vec2d center;
			unsigned int n_segments;
			std::vector<Vec2d> segments_lla;
			std::vector<Vec2d> segments;
			double area;
			MapGrid map_grid;
		};

		struct ship{
			ship(){}
			poly polygon;
			double heading;
			double length;
			double width;
		};
		
		GeometricConstr();
		GeometricConstr(int nbrCShips);
		GeometricConstr(int nbrCShips, int nbrAShips);
		
		void loadPolygon(std::string filename, vector<poly> *polygons, unsigned int *n_polygons);
		void loadHarbPolygon(std::string filename);
		void loadObstPolygon(std::string filename);
		void loadSpeePolygon(std::string filename);
		void loadGuidPolygon(std::string filename);
		void loadHarbPolygon(std::string filename, std::string gridfilename);
		void loadObstPolygon(std::string filename, std::string gridfilename);
		
		void initializeShip(ship *s, double length, double width);
		void initializeCShip(int i, double length, double width);
		void initializeAShip(int i, double length, double width);
		
		void updateShip(ship *s, double center_x, double center_y, double heading);
		void updateCShip(int i, double center_x, double center_y, double heading);
		void updateAShip(int i, double center_x, double center_y, double heading);
		
		void DistArea( poly *p1, poly *p2,                    double* g1, double* g2 );
		void DistArea( poly *p1, poly *p2, MapGrid *map_grid, double* g1, double* g2);
		void ConstrCShipCShip( unsigned int k, unsigned int l, double* g1, double* g2, bool withGrid );
		void ConstrCShipAShip( unsigned int k, unsigned int l, double* g1, double* g2, bool withGrid );
		void ConstrCShipHarbr( unsigned int k, unsigned int l, double* g1, double* g2, bool withGrid );
		void ConstrCShipObsta( unsigned int k, unsigned int l, double* g1, double* g2, bool withGrid );
		
		double SpeedBound( unsigned int k, vector<double> speedbounds );
		
		void sqDistShipGuide( unsigned int k, double* g1 );
		std::vector<Vec2d>* getHarborSegements(unsigned int l);
		std::vector<Vec2d>* getObstaclesSegements(unsigned int l);
		unsigned int getNbrObstacles();
		double getAFactor();
		double getAreaOfShip(unsigned int i);

		double lat_max;
		double lat_min;
		double lon_max;
		double lon_min;
		dsl::TypedMap* occ_map;
		void initializeOccMap(Eigen::Vector2d nr_cells, double buffer);
 		bool create_new_occ_map;
		string fOccMap;
		string tOccMap;
		Eigen::Vector2d occ_map_cell_dim_ned;
		
// 	private:
		double lat0;
		double lon0;
		double hgt0;
		
		const double EPSILON = 1.0E-5;
		const double n_ship_segments = 4;
		const double A_factor=1.0E-4;//factor to adjust A^2 to dist^2
		
		unsigned int nbrCShips;
		unsigned int nbrAShips;
		

		std::vector<ship> Cships;
		std::vector<ship> Aships;
		
		unsigned int    n_harb_polygons;
		std::vector<poly> harb_polygons;
		
		unsigned int    n_obst_polygons;
		std::vector<poly> obst_polygons;
		
		unsigned int    n_spee_polygons;
		std::vector<poly> spee_polygons;
		
		unsigned int    n_guid_polygons;
		std::vector<poly> guid_polygons;
		

		void	transformPolygonsToNED(vector<poly> *polygons, unsigned int *n_polygons);
		void 	transformHarbPolygonsToNED();
		void	transformObstPolygonsToNED();
		void 	transformGuidPolygonsToNED();
		void 	transformSpeePolygonsToNED();
		
		double 	areaOfPoly(std::vector<Vec2d>& segs);
		bool 	isPointInsidePolygon( std::vector<Vec2d>& segs, double pt_x, double pt_y );
		double 	sqDinstanceSegmentSegment( Vec2d &p1, Vec2d &q1, Vec2d &p2, Vec2d &q2 );
		double 	Clamp( double n, double min, double max );
		int 	left_of(Vec2d& a, Vec2d& b, Vec2d& c);
		int 	line_sect(Vec2d& x0, Vec2d& x1, Vec2d& y0, Vec2d& y1, Vec2d& res);
		int 	poly_winding(std::vector<Vec2d>& p);
		void 	poly_edge_clip(std::vector<Vec2d>& sub, Vec2d& x0, Vec2d& x1, int left, std::vector<Vec2d>& res);
		void	poly_clip(std::vector<Vec2d>& sub, std::vector<Vec2d>& clip, std::vector<Vec2d>& p2);
		
		double 	sqDinstancePolyPoly( std::vector<Vec2d>& segs1, std::vector<Vec2d>& segs2 );
		double 	sqDinstancePolySegment( std::vector<Vec2d>& segs1, std::vector<Vec2d>& segs2, unsigned int segs2_idx );
		double 	sqDistancePolyPoint( std::vector<Vec2d>& segs1, Vec2d &p2 );
		double 	sqDinstanceSegmentPoint( Vec2d &p1, Vec2d &q1, Vec2d &p2 );

};

#endif
