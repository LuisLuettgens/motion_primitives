#include "dp_0xA1.h"

using namespace std;


dp_0xA1::dp_0xA1 ( ) {
	NBR_ENTRIES_IN_HEADER  =  8; 
	NBR_ENTRIES_IN_PCKT    =  4; 
	NBR_ENTRIES_IN_FOOTER  =  1; 
	NBR_OF_PCKTS           =  1;
	
	header .resize(NBR_ENTRIES_IN_HEADER);
	packInt.resize(NBR_OF_PCKTS);
	pack   .resize(NBR_OF_PCKTS);
	footer .resize(NBR_ENTRIES_IN_FOOTER);

	header[0] = {"Preamble"          , "Preamble"  , 0};
	header[1] = {"Packet Identifier" , "Packet ID" , 0};
	header[2] = {"Product ID"        , "Ship ID"   , 0};
	header[3] = {"Notaus ID"         , "Notaus ID" , 0};
	

	for(int i=0; i<NBR_OF_PCKTS; i++) {
		packInt[i].resize(NBR_ENTRIES_IN_PCKT);
		pack[i]   .resize(NBR_ENTRIES_IN_PCKT);

		packInt[i][ 0] = {"Position Latitude"         , "Latitude"   , 0};
		packInt[i][ 1] = {"Position Longitude"        , "Longitude"  , 0};
		packInt[i][ 2] = {"Statusbyte"                , "Statusbyte" , 0};
		packInt[i][ 3] = {"Time"                      , "Time"       , 0};
	}
	
	footer[0] = {"RC-16 Checksum", "Checksum", 0};
	
	SIZE_OF_MSG = get_size_of_msg();
	msg =  new char[ SIZE_OF_MSG ];
}

void dp_0xA1::toTW(int k, double lat0, double lon0, double h0) {
}
