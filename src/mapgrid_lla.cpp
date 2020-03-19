#include "mapgrid_lla.h"


using namespace std;


MapGrid::MapGrid(std::string xml_filename, double lat0, double lon0, double h0):active(0)
{
	read_xml(xml_filename, lat0, lon0, h0);
}

void MapGrid::read_xml(std::string xml_filename, double lat0, double lon0, double h0)
{

	XMLParser parser;
	XMLNode *node = parser.Parse(xml_filename);
	if(node){
		XMLNode *grid = node->GetFirstChild("grid");
		if (grid) {
			
			XMLNode *n = grid->GetFirstChild("width");
			if (n) {
				width = atof(n->GetText().c_str());
			}else{
				std::cout<<"Error in MapGrid::read_xml: cannot open node: width"<<std::endl;
				exit(1);
			}
			
			n = grid->GetFirstChild("lat");
			if (n) {
				std::string str = n->GetText();
				std::string val= "";
				unsigned int i=1;
				while ( i<str.size() && str.at ( i ) !=']' ) {
					while( i<str.size() && str.at ( i ) !=' ' ){
						val += str.at ( i );
						i++;
					}
					lat.push_back( atof(val.c_str()) );
					val= "";
					i++;
				}
			}else{
				std::cout<<"Error in MapGrid::read_xml: cannot open node: lat"<<std::endl;
				exit(1);
			}
			
			n = grid->GetFirstChild("lon");
			if (n) {
				std::string str = n->GetText();
				std::string val= "";
				unsigned int i=1;
				while ( i<str.size() && str.at ( i ) !=']' ) {
					while( i<str.size() && str.at ( i ) !=' ' ){
						val += str.at ( i );
						i++;
					}
					lon.push_back( atof(val.c_str()) );
					val= "";
					i++;
				}
			}else{
				std::cout<<"Error in MapGrid::read_xml: cannot open node: lon"<<std::endl;
				exit(1);
			}
			
			n = grid->GetFirstChild("min_lat");
			if (n) {
				min_lat = atof(n->GetText().c_str());
			}else{
				std::cout<<"Error in MapGrid::read_xml: cannot open node: min_lat"<<std::endl;
				exit(1);
			}
			
			n = grid->GetFirstChild("max_lat");
			if (n) {
				max_lat = atof(n->GetText().c_str());
			}else{
				std::cout<<"Error in MapGrid::read_xml: cannot open node: max_lat"<<std::endl;
				exit(1);
			}
			
			n = grid->GetFirstChild("min_lon");
			if (n) {
				min_lon = atof(n->GetText().c_str());
			}else{
				std::cout<<"Error in MapGrid::read_xml: cannot open node: min_lon"<<std::endl;
				exit(1);
			}
			
			n = grid->GetFirstChild("max_lon");
			if (n) {
				max_lon = atof(n->GetText().c_str());
			}else{
				std::cout<<"Error in MapGrid::read_xml: cannot open node: max_lon"<<std::endl;
				exit(1);
			}
			
			n = grid->GetFirstChild("lat_size");
			if (n) {
				lat_size = atof(n->GetText().c_str());
			}else{
				std::cout<<"Error in MapGrid::read_xml: cannot open node: lat_size"<<std::endl;
				exit(1);
			}
			
			n = grid->GetFirstChild("lon_size");
			if (n) {
				lon_size = atof(n->GetText().c_str());
			}else{
				std::cout<<"Error in MapGrid::read_xml: cannot open node: lon_size"<<std::endl;
				exit(1);
			}
			
			harbor_area.resize(lat_size);
			for (unsigned int i=0; i<lat_size; i++) {
				harbor_area[i].resize(lon_size);
			}
			unsigned int count_x=0;
			unsigned int count_y=0;
			n = grid->GetFirstChild("harbor_area");
			if (n) {
				XMLNode *m = n->GetFirstChild("tr");
				while (m) {
					XMLNode *l = m->GetFirstChild("td");
					while (l) {
						XMLNode *kx = l->GetFirstChild("lat");
						XMLNode *ky = l->GetFirstChild("lon");
						if (kx && ky) {
							
							std::string strx = kx->GetText();
							std::string stry = ky->GetText();
							std::string valx= "";
							std::string valy= "";
							unsigned int i=1;
							unsigned int j=1;
							while ( i<strx.size() && strx.at ( i ) !=']'
									&& j<stry.size() && stry.at ( j ) !=']'
							) {
								while( i<strx.size() && strx.at ( i ) !=' ' ){
									valx += strx.at ( i );
									i++;
								}
								while( j<stry.size() && stry.at ( j ) !=' ' ){
									valy += stry.at ( j );
									j++;
								}
								
								double x;
								double y;
								double z;
								LLAtoNED(atof(valx.c_str()), atof(valy.c_str()), 0, lat0, lon0, h0, &x, &y, &z);
								harbor_area[count_x][count_y].push_back(Vec2d(x, y ) );
								valx= "";
								valy= "";
								i++;
								j++;
								
							}
							
						}
						l = m->GetNextChild("td");
						count_x++;
					}
					m = n->GetNextChild("tr");
					count_x=0;
					count_y++;
				}
			}else{
				std::cout<<"Error in MapGrid::read_xml: cannot open node: harbor_area"<<std::endl;
				exit(1);
			}
			
			segs_in_dist.resize(lat_size);
			for (unsigned int i=0; i<lat_size; i++) {
				segs_in_dist[i].resize(lon_size);
			}
			count_x=0;
			count_y=0;
			n = grid->GetFirstChild("segs_in_dist");
			if (n) {
				XMLNode *m = n->GetFirstChild("tr");
				while (m) {
					XMLNode *l = m->GetFirstChild("td");
					while (l) {
						std::string str = l->GetText();
						if ( !str.empty() ) {
							std::string val= "";
							unsigned int i=1;
							while ( i<str.size() && str.at ( i ) !=']') {
								while( i<str.size() && str.at ( i ) !=' ' ){
									val += str.at ( i );
									i++;
								}
								segs_in_dist[count_x][count_y].push_back(atof( val.c_str())-1 );//-1: converting from matlab indexing to c++
								val= "";
								i++;
							}
						}
						l = m->GetNextChild("td");
						count_x++;
					}
					m = n->GetNextChild("tr");
					count_x=0;
					count_y++;
				}
			}else{
				std::cout<<"Error in MapGrid::read_xml: cannot open node: segs_in_dist"<<std::endl;
				exit(1);
			}
			
		}else{
			std::cout<<"Error in MapGrid::read_xml: cannot open node: grid"<<std::endl;
			exit(1);
		}
	}else{
		std::cout<<"Error in MapGrid::read_xml: cannot open xml-file: "<<xml_filename<<std::endl;
		exit(1);
	}
	active = 1;
}



