#include "dp_0xA3.h"

using namespace std;


dp_0xA3::dp_0xA3 ( ) {
	
	NBR_ENTRIES_IN_HEADER  =   6; 
	NBR_ENTRIES_IN_PCKT    =  14; 
	NBR_ENTRIES_IN_FOOTER  =   1; 
	NBR_OF_PCKTS           = 600;
	
	header .resize(NBR_ENTRIES_IN_HEADER);
	packInt.resize(NBR_OF_PCKTS);
	pack   .resize(NBR_OF_PCKTS);
	footer .resize(NBR_ENTRIES_IN_FOOTER);

	header[0] = {"Preamble"             , "Preamble"     , 0};
	header[1] = {"Packet Identifier"    , "Packet ID"    , 0};
	header[2] = {"Product ID"           , "Product ID"   , 0};
	header[3] = {"Trajectory ID"        , "Trajectory ID", 0};
	header[4] = {"Sub-packet Identifier", "Sub-packet I" , 0};
	header[5] = {"Counter, Heartbeat"   , "Counter"      , 0};

	for(int i=0; i<NBR_OF_PCKTS; i++) {
		packInt[i].resize(NBR_ENTRIES_IN_PCKT);
		pack[i]   .resize(NBR_ENTRIES_IN_PCKT);
		
		packInt[i][ 0] = {"Time"                     , "Time(s)"       , 0};
		packInt[i][ 1] = {"Position Latitude"        , "Lat"           , 0};
		packInt[i][ 2] = {"Position Longitude"       , "Lon"           , 0};
		packInt[i][ 3] = {"Absolute angle"           , "Hdg(deg)"      , 0};
		packInt[i][ 4] = {"Translational velocity"   , "u-vel.(m/s)"   , 0};
		packInt[i][ 5] = {"Translational velocity"   , "v-vel.(m/s)"   , 0};
		packInt[i][ 6] = {"Rotational velocity"      , "RoT(deg/s)"    , 0};
		packInt[i][ 7] = {"Rudder deflection"        , "Rudder1(deg)"  , 0};
		packInt[i][ 8] = {"Rudder deflection"        , "Rudder2(deg)"  , 0};
		packInt[i][ 9] = {"Revolutions of Propulsion", "Prop.1(rpm)"   , 0};
		packInt[i][10] = {"Revolutions of Propulsion", "Prop.2(rpm)"   , 0};
		packInt[i][11] = {"Revolutions of Propulsion", "Thruster1(rpm)", 0};
		packInt[i][12] = {"Revolutions of Propulsion", "Thruster2(rpm)", 0};
		packInt[i][13] = {"Control Configuration"    , "Ctrl Config"   , 0};
	}
	
	footer[0] = {"RC-16 Checksum", "Checksum", 0};
	
	SIZE_OF_MSG = get_size_of_msg();
	msg =  new char[ SIZE_OF_MSG ];
}

