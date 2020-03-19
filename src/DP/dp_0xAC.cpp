#include "dp_0xAC.h"

using namespace std;


dp_0xAC::dp_0xAC ( ) {
	
	NBR_ENTRIES_IN_HEADER  =   5; 
	NBR_ENTRIES_IN_PCKT    =  29; 
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
	header[4] = {"Nanoseconds"       , "Nanosec"   , 0};

	for(int i=0; i<NBR_OF_PCKTS; i++) {
		packInt[i].resize(NBR_ENTRIES_IN_PCKT);
		pack[i]   .resize(NBR_ENTRIES_IN_PCKT);

		packInt[i][ 0] = {"Product ID"               , "Ship ID"       , 0};
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
		
		packInt[i][14] = {"Absolute angle", "TWD"  , 0};
		packInt[i][15] = {"Absolute speed", "TWS"  , 0};
		packInt[i][16] = {"Absolute angle", "Set"  , 0};
		packInt[i][17] = {"Absolute speed", "Drift", 0};
		
		packInt[i][18] = {"Integrity Flag", "Integrity of Data", 0};
		
		packInt[i][19] = {"Position Latitude"     , "Lat-Std"        , 0};
		packInt[i][20] = {"Position Longitude"    , "Lon-Std"        , 0};
		packInt[i][21] = {"Absolute angle"        , "Hdg(deg)-Std"   , 0};
		packInt[i][22] = {"Translational velocity", "u-vel.(m/s)-Std", 0};
		packInt[i][23] = {"Translational velocity", "v-vel.(m/s)-Std", 0};
		packInt[i][24] = {"Rotational velocity"   , "RoT(deg/s)-Std" , 0};
		
		packInt[i][25] = {"Absolute angle", "TWD-Std"  , 0};
		packInt[i][26] = {"Absolute speed", "TWS-Std"  , 0};
		packInt[i][27] = {"Absolute angle", "Set-Std"  , 0};
		packInt[i][28] = {"Absolute speed", "Drift-Std", 0};
	}
	
	footer[0] = {"RC-16 Checksum", "Checksum", 0};
	
	SIZE_OF_MSG = get_size_of_msg();
	msg =  new char[ SIZE_OF_MSG ];
}

void dp_0xAC::toTW(double lat0, double lon0, double h0) {
	
	for(int k=0; k<ph0->nbrCShips; k++) {
		int idx = prodIDmap[packInt[k][0].value];
		
		double lat   = pack[k][1];
		double lon   = pack[k][2];
		double psi   = pack[k][3]*M_PI/180;
		double u     = pack[k][4];
		double v     = pack[k][5];
		double r     = pack[k][6]*M_PI/180;
		double delta = pack[k][7]*M_PI/180;
		double n     = pack[k][9]/60;
		
		cout << std::setprecision(10);
		for(int i=1; i<10; i++) cout << pack[k][i] << ", ";
		cout << endl;
		
		double x0, y0, z0;
		LLAtoNED(lat, lon, 0, lat0, lon0, h0, &x0, &y0, &z0);
		
		
		ph0->var_b_init[idx][0] = x0;
		ph0->var_b_init[idx][1] = y0;
		ph0->var_b_init[idx][2] = psi;
		ph0->var_b_init[idx][3] = u;
		ph0->var_b_init[idx][4] = v;
		ph0->var_b_init[idx][5] = r;
		ph0->var_b_init[idx][6] = delta;
		ph0->var_b_init[idx][7] = n;
		
		for(int i=0; i<8; i++) {
			ph0->var_b_init[idx][i] /= ph0->scale[idx][i];
		}
	}
}

