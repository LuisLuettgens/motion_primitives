#include "dp_0xAD.h"

using namespace std;


dp_0xAD::dp_0xAD():dp_base() {

	NBR_ENTRIES_IN_HEADER  =   4; 
	NBR_ENTRIES_IN_PCKT    =  11; 
	NBR_ENTRIES_IN_FOOTER  =   1; 
	NBR_OF_PCKTS           =  10;
	
	header .resize(NBR_ENTRIES_IN_HEADER);
	packInt.resize(NBR_OF_PCKTS);
	pack   .resize(NBR_OF_PCKTS);
	footer .resize(NBR_ENTRIES_IN_FOOTER);
	
	header[0] = {"Preamble"          , "Preamble"  , 0};
	header[1] = {"Packet Identifier" , "Packet ID" , 0};
	header[2] = {"Number"            , "Nbr Ships" , 0};
	header[3] = {"Time"              , "Time"      , 0};

	for(int i=0; i<NBR_OF_PCKTS; i++) {
		packInt[i].resize(NBR_ENTRIES_IN_PCKT);
		pack[i]   .resize(NBR_ENTRIES_IN_PCKT);

		packInt[i][ 0] = {"Product ID"          , "Ship ID"     , 0};
		packInt[i][ 1] = {"Euclidean distance"  , "Length"      , 0};
		packInt[i][ 2] = {"Euclidean distance"  , "Width"       , 0};
		packInt[i][ 3] = {"Position Latitude"   , "Lat"         , 0};
		packInt[i][ 4] = {"Position Longitude"  , "Lon"         , 0};
		packInt[i][ 5] = {"Absolute angle"      , "Hdg(deg)"    , 0};
		packInt[i][ 6] = {"Absolute angle"      , "COG(m/s)"    , 0};
		packInt[i][ 7] = {"Absolute speed"      , "SOG(m/s)"    , 0};
		packInt[i][ 8] = {"Rotational velocity" , "RoT(deg/s)"  , 0};
		packInt[i][ 9] = {"Distance (short)"    , "Position Acc", 0};
		packInt[i][10] = {"Integrity Flag"      , "Integrity"   , 0};
	}
	
	footer[0] = {"RC-16 Checksum", "Checksum", 0};
	
	SIZE_OF_MSG = get_size_of_msg();
	msg =  new char[ SIZE_OF_MSG ];
}

void dp_0xAD::toTW(double lat0, double lon0, double h0) {
	
	ph0->nbrAShips = header[2].value;
	ph0->AISShips.resize(ph0->nbrAShips);
	
	ph0->geom_constr.nbrAShips = ph0->nbrAShips;
	ph0->geom_constr.Aships.resize(ph0->nbrAShips);
	
	for(int k=0; k<ph0->nbrAShips; k++) {
		ph0->geom_constr.Aships[k].polygon.n_segments = ph0->geom_constr.n_ship_segments;
		ph0->geom_constr.Aships[k].polygon.segments_lla.resize(ph0->geom_constr.n_ship_segments);
		ph0->geom_constr.Aships[k].polygon.segments    .resize(ph0->geom_constr.n_ship_segments);
		ph0->geom_constr.initializeAShip(k, pack[k][1], pack[k][2]);

		double lat   = pack[k][ 3];
		double lon   = pack[k][ 4];
		double psi   = pack[k][ 5]*M_PI/180;
		double cog   = pack[k][ 6]*M_PI/180;
		double sog   = pack[k][ 7];
		double r     = pack[k][ 8]*M_PI/180;
		double acc   = pack[k][ 9];
		double integ = pack[k][10];
		
		double x0, y0, z0;
		LLAtoNED(lat, lon, 0, lat0, lon0, h0, &x0, &y0, &z0);
		
		ph0->AISShips[k].x0        = x0;
		ph0->AISShips[k].y0        = y0;
		ph0->AISShips[k].psi       = psi;
		ph0->AISShips[k].cog       = cog;
		ph0->AISShips[k].sog       = sog;
		ph0->AISShips[k].r         = r;
		ph0->AISShips[k].accuracy  = acc;
		ph0->AISShips[k].integrity = integ;
	}
}

