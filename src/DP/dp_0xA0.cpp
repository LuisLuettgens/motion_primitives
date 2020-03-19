#include "dp_0xA0.h"

using namespace std;


dp_0xA0::dp_0xA0 ( ) {
	NBR_ENTRIES_IN_HEADER  =  7; 
	NBR_ENTRIES_IN_PCKT    =  1; 
	NBR_ENTRIES_IN_FOOTER  =  1; 
	NBR_OF_PCKTS           =  1;
	
	header .resize(NBR_ENTRIES_IN_HEADER);
	packInt.resize(NBR_OF_PCKTS);
	pack   .resize(NBR_OF_PCKTS);
	footer .resize(NBR_ENTRIES_IN_FOOTER);

	header[0] = {"Preamble"          , "Preamble"      , 0};
	header[1] = {"Packet Identifier" , "Packet ID"     , 0};
	header[2] = {"Counter, Heartbeat", "Counter"       , 0};
	header[3] = {"Product ID"        , "Product ID"    , 0};
	header[4] = {"Product ID"        , "Product ID"    , 0};
	header[5] = {"Alert ID"          , "Alert ID"      , 0};
	header[6] = {"Alert Severity"    , "Alert Severity", 0};
	

	for(int i=0; i<NBR_OF_PCKTS; i++) {
		packInt[i].resize(NBR_ENTRIES_IN_PCKT);
		pack[i]   .resize(NBR_ENTRIES_IN_PCKT);

		packInt[i][0] = {"Time", "Time", 0};
	}
	
	footer[0] = {"RC-16 Checksum", "Checksum", 0};
	
	SIZE_OF_MSG = get_size_of_msg();
	msg =  new char[ SIZE_OF_MSG ];
}

void dp_0xA0::toTW(int k, double lat0, double lon0, double h0) {
}

