#ifndef mapgrid_lla_h
#define mapgrid_lla_h

#include "utils/functions.h"


class MapGrid {
	public:
		MapGrid():active(0){}
		MapGrid(std::string xml_filename, double lat0, double lon0, double h0);
		double width;
		std::vector<double> lat;
		std::vector<double> lon;
		unsigned int lat_size;
		unsigned int lon_size;
		double min_lat;
		double max_lat;
		double min_lon;
		double max_lon;
		std::vector< std::vector < std::vector<Vec2d> > > harbor_area;
		std::vector< std::vector < std::vector<unsigned int> > > segs_in_dist;
		bool active;
		void read_xml(std::string xml_filename, double lat0, double lon0, double h0);
};

#endif
