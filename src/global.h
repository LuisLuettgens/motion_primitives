#ifndef global_H
#define global_H

#include "optimal_traj_Nships.h"

extern OptTraj *ph0;
extern tw::TWfolder *folder0;
extern int *idx_viewer0;

extern bool first_time_in_0xAC;


extern string ipad;
extern int port1;
extern int port2;
extern size_t MaxsizeofMsg;
extern std::map<string, struct Range> DTranges;
extern std::map<string, struct Description> DD;
extern std::map<int, struct AISShipDim> AISShipDims;
extern std::map<uint16_t,int> prodIDmap;

#endif

