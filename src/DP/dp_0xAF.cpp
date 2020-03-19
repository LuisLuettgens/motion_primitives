#include "dp_0xAF.h"

using namespace std;


dp_0xAF::dp_0xAF ( ) {
	
	NBR_ENTRIES_IN_HEADER  =  3; 
	NBR_ENTRIES_IN_PCKT    =  7; 
	NBR_ENTRIES_IN_FOOTER  =  1; 
	NBR_OF_PCKTS           =  1;

	header .resize(NBR_ENTRIES_IN_HEADER);
	packInt.resize(NBR_OF_PCKTS);
	pack   .resize(NBR_OF_PCKTS);
	footer .resize(NBR_ENTRIES_IN_FOOTER);

	header[0] = {"Preamble"          , "Preamble"  , 0};
	header[1] = {"Packet Identifier" , "Packet ID" , 0};
	header[2] = {"Product ID"        , "Ship ID"   , 0};

	for(int i=0; i<NBR_OF_PCKTS; i++) {
		packInt[i].resize(NBR_ENTRIES_IN_PCKT);
		pack[i]   .resize(NBR_ENTRIES_IN_PCKT);

		packInt[i][ 0] = {"Position Latitude"         , "Lat"            , 0};
		packInt[i][ 1] = {"Position Longitude"        , "Lon"            , 0};
		packInt[i][ 2] = {"Absolute angle"            , "Hdg(deg)"       , 0};
		packInt[i][ 3] = {"Translational velocity"    , "u-vel.(m/s)"    , 0};
		packInt[i][ 4] = {"Translational velocity"    , "v-vel.(m/s)"    , 0};
		packInt[i][ 5] = {"Rotational velocity"       , "RoT(deg/s)"     , 0};
		packInt[i][ 6] = {"Final States Specificatons", "States Spec", 0};
	}
	
	footer[0] = {"RC-16 Checksum", "Checksum", 0};
	
	SIZE_OF_MSG = get_size_of_msg();
	msg =  new char[ SIZE_OF_MSG ];
}

void dp_0xAF::toTW(int k, double lat0, double lon0, double h0) {

	double lat   = pack[0][0];
	double lon   = pack[0][1];
	double psi   = pack[0][2]*M_PI/180;
	double u     = pack[0][3];
	double v     = pack[0][4];
	double r     = pack[0][5]*M_PI/180;
	
	int nbr = packInt[0][6].value;
	bool bin_dig[8];
	for(int i=0; i<8; i++) {
		bin_dig[7-i] = nbr/(int)pow(2,7-i); 
		nbr = nbr%(int)pow(2,7-i);
	}
	
	double x0, y0, z0;
	LLAtoNED(lat, lon, 0, lat0, lon0, h0, &x0, &y0, &z0);
	
	ph0->final_s[k][0] = bin_dig[0] ? x0  : nan(""); 
	ph0->final_s[k][1] = bin_dig[1] ? y0  : nan(""); 
	ph0->final_s[k][2] = bin_dig[2] ? psi : nan(""); 
	ph0->final_s[k][3] = bin_dig[3] ? u   : nan(""); 
	ph0->final_s[k][4] = bin_dig[4] ? v   : nan(""); 
	ph0->final_s[k][5] = bin_dig[5] ? r   : nan(""); 
	ph0->final_s[k][6] = nan(""); 
	ph0->final_s[k][7] = nan(""); 
	ph0->final_s[k][8] = nan(""); 

	for(int i=0; i<9; i++) ph0->final_s[k][i] /= ph0->scale[k][i];
}

