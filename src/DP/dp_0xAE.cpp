#include "dp_0xAE.h"

using namespace std;


dp_0xAE::dp_0xAE ( ) {
	
	NBR_ENTRIES_IN_HEADER  =   4; 
	NBR_ENTRIES_IN_PCKT    =   6; 
	NBR_ENTRIES_IN_FOOTER  =   1; 
	NBR_OF_PCKTS           =  10;
	
	header .resize(NBR_ENTRIES_IN_HEADER);
	packInt.resize(NBR_OF_PCKTS);
	pack   .resize(NBR_OF_PCKTS);
	footer .resize(NBR_ENTRIES_IN_FOOTER);

	header[0] = {"Preamble"          , "Preamble"      , 0};
	header[1] = {"Packet Identifier" , "Packet ID"     , 0};
	header[2] = {"Number"            , "Nbr Obstacles" , 0};
	header[3] = {"Time"              , "Time"          , 0};

	for(int i=0; i<NBR_OF_PCKTS; i++) {
		packInt[i].resize(NBR_ENTRIES_IN_PCKT);
		pack[i]   .resize(NBR_ENTRIES_IN_PCKT);

		packInt[i][ 0] = {"Position Latitude" , "Lat"       , 0};
		packInt[i][ 1] = {"Position Longitude", "Lon"       , 0};
		packInt[i][ 2] = {"Euclidean distance", "Size"      , 0};
		packInt[i][ 3] = {"Position Latitude" , "Lat-Std" , 0};
		packInt[i][ 4] = {"Position Longitude", "Lon-Std" , 0};
		packInt[i][ 5] = {"Euclidean distance", "Size-Std", 0};
	}
	
	footer[0] = {"RC-16 Checksum", "Checksum", 0};
	
	SIZE_OF_MSG = get_size_of_msg();
	msg =  new char[ SIZE_OF_MSG ];
}

void dp_0xAE::toTW() {
	
	ph0->geom_constr.obst_polygons.resize(ph0->nbrObstacles_init);
	ph0->geom_constr.n_obst_polygons = ph0->nbrObstacles_init;
	
	for(int k=0; k<header[2].value; k++) {
		ofstream outfile;
		outfile.open("obs_temp.txt");
		outfile << pack[k][1] << " " << pack[k][0] << endl;
		outfile.close();
		
		ph0->geom_constr.loadObstPolygon("obs_temp.txt", "");
		ph0->nbrObstacles = ph0->geom_constr.getNbrObstacles();
		ph0->geom_constr.transformObstPolygonsToNED();
	}
}

