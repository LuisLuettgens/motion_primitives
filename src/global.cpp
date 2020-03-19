#include "global.h"

#include "optimal_traj_Nships.h"

OptTraj  *ph0    =0;
tw::TWfolder *folder0=0;
int *idx_viewer0 = 0;


string ipad = "127.0.0.1";
int port1 = 1502;
int port2 = 1503;

bool first_time_in_0xAC = true;

size_t MaxsizeofMsg = 20000;

std::map<string, struct Range> DTranges = {
	{ "uInt8"  , {1,           0,        255}},
	{ "Int8"   , {1,        -128,        127}},
	{ "uInt16" , {2,           0,      65535}},
	{ "Int16"  , {2,      -32768,      32767}},
	{ "uInt32" , {4,           0, 4294967295}},
	{ "Int32"  , {4, -2147483648, 2147483647}},
};

std::map<string, struct Description> DD = {
	{ "Position Latitude"         , {1  , "Int32" ,   -180,        180}},
	{ "Position Longitude"        , {2  , "Int32" ,   -180,        180}},
	{ "Distance (short)"          , {3  , "Int16" ,  -1000,       1000}},
	{ "Euclidean distance"        , {5  , "uInt16" ,     0,     655.35}},
	{ "Absolute angle"            , {7  , "uInt16",      0,        360}},
	{ "Absolute speed"            , {8  , "uInt16",      0,        256}},
	{ "Translational velocity"    , {9  , "Int16" ,   -256,        256}},
	{ "Rotational velocity"       , {10 , "Int16" ,   -650,        650}},
	{ "Revolutions of Propulsion" , {31 , "Int16" , -10000,      10000}},
	{ "Rudder deflection"         , {51 , "Int16" ,    -45,         45}},
	{ "Product ID"                , {150, "uInt16",      0,      65535}},
	{ "Alert ID"                  , {151, "uInt16",      0,      65535}},
	{ "Alert Severity"            , {152, "uInt8" ,      0,        255}},
	{ "Notaus ID"                 , {153, "uInt8" ,      0,        255}},
	{ "Trajectory ID"             , {154, "uInt16",      0,      65535}},
	{ "Sub-packet Identifier"     , {155, "uInt8" ,      1,        256}},
	{ "Control Configuration"     , {237, "uInt8" ,      0,        255}},
	{ "Boolean Flag"              , {238, "uInt8" ,      0,        255}},
	{ "Integrity Flag"            , {239, "uInt8" ,      0,        255}},
	{ "Final States Specificatons", {241, "uInt8" ,      0,        255}},
	{ "Number"                    , {247, "uInt8" ,      0,        255}},
	{ "Time"                      , {248, "uInt32",      0, 4294967295}},
	{ "Nanoseconds"               , {249, "uInt32",      0,  999999999}},
	{ "Statusbyte"                , {251, "uInt8" ,      0,        255}},
	{ "Preamble"                  , {253, "uInt16",      0,          0}},
	{ "Counter, Heartbeat"        , {254, "uInt32",      0, 4294967295}},
	{ "RC-16 Checksum"            , {255, "uInt16",      0,          0}},
	{ "Packet Identifier"         , {256, "uInt8" ,      0,          0}},
	{ "Integrity Flag"            , {257, "uInt8" ,      0,        255}},
};

std::map<int, struct AISShipDim> AISShipDims = {
	{ 0  , {200., 47.5}},
	{ 1  , {200., 47.5}},
};

std::map<uint16_t,int> prodIDmap = {
	{0, 0}, 
	{1, 1}, 
	{2, 2}
};



