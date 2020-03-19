/*    Optimal Trajectory Calculations     */

#ifdef WIN32
#include "windows.h"
#endif
#include "optimal_traj_Nships.h"
#include "gridpathsearch.h"
#include "DP/dp_0xA0.h"
#include "DP/dp_0xA1.h"
#include "DP/dp_0xA2.h"
#include "DP/dp_0xA3.h"
#include "DP/dp_0xA4.h"
#include "DP/dp_0xAC.h"
#include "DP/dp_0xAD.h"
#include "DP/dp_0xAE.h"
#include "DP/dp_0xAF.h"

#include <fstream>
#include <time.h>

 
void log( const std::string &text )
{
//     std::ofstream log_file(
//         "log_file.txt", std::ios_base::out | std::ios_base::app );
//     log_file << text << endl;
}

bool is_all_true(vector<bool> bool_vec) {
	int sum=0;
	for(int k=0; k<bool_vec.size(); k++) {
		sum += bool_vec[k];
	}
	if(sum==bool_vec.size()) return true;
	return false;
}

bool is_all_false(vector<bool> bool_vec) {
	int sum=0;
	for(int k=0; k<bool_vec.size(); k++) {
		sum += bool_vec[k];
	}
	if(sum==0) return true;
	return false;
}

TypeOfVehicle vehicle = SHIP;

int main ( int argv, char* argc[] ) {
	
	string transxml = "transworhp.xml";
	string shipvxml;
	
	vector<string> harb_f;
	vector<string> grid_f;
	vector<string> guid_f;
	vector<string> obst_f;
	vector<string> obst_grid_f;
	vector<string> spee_f;
	
	double lat0;
	double lon0;
	double hgt0;
	
	string fOccMap;
	string tOccMap;
	
	vector <string> coop_ships_IDs;
	
	int nbrAShips;
	int nbrAShipsMax;
	int nbrObstacles;
	int nbrObstacles_init;
	int nbrExtraObstacles;
	int nbrObstaclesMax;
	
	vector<double>         distSH;
	vector<double>         distSO;
	vector<vector<double>> distSA;
	vector<vector<double>> distSS;
	vector<double> speedbounds;
	
	int Ninterp;
	
	vector <double> ratioToFirst;
	
	bool withCollReg;
	EndTimeBool ENDTIME;
	vector<double> tEnd_nu ;
	bool collavoidSH ;
	bool collavoidSO;
	bool collavoidSA;
	bool collavoidSS;
	bool speedLimit;
	
	vector<double>         tWgt;
	vector<double>         eWgt;
	vector<vector<double>> fWgt;
	vector<vector<double>> cWgt;
	vector<double>         plow_nu;
	vector<double>         pupp_nu;
	
	vector<double>          scale_t;
	vector <vector<double>> scale;
	vector <vector<double>> scale_ctrl;

	double viewer_t_scale;
	
	double stepSize;
	
	double buffer;
	vector<double> ned_map_delta;
	vector<double> ned_grid_delta;
	
	bool soft_final_constr;
	
	if(vehicle == SHIP) {
		shipvxml = "shipviewer.xml";
		
		lat0 = 54.17057475;
		lon0 = 12.10074142;
		hgt0 = 0;
		
		harb_f.push_back("RostockHarbor/harbor_simplified.txt");
		grid_f.push_back("RostockHarbor/harbor_simplified.xml");
		obst_f.push_back("RostockHarbor/obstacle_01_simplified.txt");
		obst_f.push_back("RostockHarbor/obstacle_02_simplified.txt");
		obst_f.push_back("RostockHarbor/obstacle_03_simplified.txt");
		obst_f.push_back("RostockHarbor/obstacle_04_simplified.txt");
		obst_f.push_back("RostockHarbor/obstacle_05_simplified.txt");
		obst_f.push_back("RostockHarbor/obstacle_06_simplified.txt");
		obst_f.push_back("RostockHarbor/obstacle_07_simplified.txt");
		obst_f.push_back("RostockHarbor/obstacle_08_simplified.txt");
		obst_f.push_back("RostockHarbor/obstacle_09_simplified.txt");
		obst_f.push_back("RostockHarbor/obstacle_10_simplified.txt");
		obst_grid_f.push_back("");
		obst_grid_f.push_back("");
		obst_grid_f.push_back("");
		obst_grid_f.push_back("");
		obst_grid_f.push_back("");
		obst_grid_f.push_back("");
		obst_grid_f.push_back("");
		obst_grid_f.push_back("");
		obst_grid_f.push_back("");
		obst_grid_f.push_back("");
		guid_f.push_back("guide.txt");
		spee_f.push_back("speedbound.txt");
		
		fOccMap = "occ_map_braunschweig.ppm";
		tOccMap = "occ_map_braunschweig.ppm";
		
		coop_ships_IDs.push_back("1");
		coop_ships_IDs.push_back("2");
		
		nbrAShips = 0;
		nbrAShipsMax = 10;
		nbrExtraObstacles = 10;
		
		distSH.push_back(50);
		distSH.push_back(50);
		
		distSO.push_back(60);
		distSO.push_back(60);
		
		distSA.push_back({100});
		distSA.push_back({100});

		distSS.push_back({100});

		speedbounds.push_back(2.0);
		
		Ninterp = 2;
		
		ratioToFirst.push_back(1.0);
		ratioToFirst.push_back(0.9);
		
		withCollReg = 0;
		ENDTIME = OPEN_WITH_FIXED_RATIOS;
		tEnd_nu.push_back(300);
		tEnd_nu.push_back(300);
		
		collavoidSH = 1;
		collavoidSO = 1;
		collavoidSA = 1;
		collavoidSS = 1;
		speedLimit  = 1;
		
		tWgt.push_back(2);
		tWgt.push_back(2);

		eWgt.push_back(2);
		eWgt.push_back(2);

		fWgt.push_back({200, 200, nan(""), 10, 10, 10, nan(""), nan("")});
		fWgt.push_back({200, 200, nan(""), 10, 10, 10, nan(""), nan("")});

		cWgt.push_back({10, 10, 0.003});
		cWgt.push_back({10, 10, 0.003});

		plow_nu.push_back(100);
		plow_nu.push_back(100);
		pupp_nu.push_back(6000);
		pupp_nu.push_back(6000);
		
		scale_t.push_back(1500);
		scale_t.push_back(700);
		
		scale.push_back({ 2000., 2000., 6.0, 3.5, 0.1, 0.001, 0.3, 1., 1. });
		scale.push_back({ 2000., 2000., 6.0, 3.5, 0.1, 0.001, 0.3, 1., 1. });

		scale_ctrl.push_back({ 0.001, 0.1 });
		scale_ctrl.push_back({ 0.02, 0.1  });
		
		viewer_t_scale = 120.0;
		
		stepSize = 250;
		
		buffer = 60.;
		ned_map_delta .push_back(20.0);
		ned_map_delta .push_back(20.0);
		ned_grid_delta.push_back(10.0);
		ned_grid_delta.push_back(10.0);
		
		soft_final_constr = 1;
	}
	else {
		shipvxml = "buggyviewer.xml";
		
		lat0 = 53.105988;
		lon0 =  8.852085;
		hgt0 =  0;
		
		harb_f.push_back("BuggyHarbor/harbor.txt");
		grid_f.push_back("");
		guid_f.push_back("");
		spee_f.push_back("");
		
		fOccMap = "occ_map_buggy__.ppm";
		tOccMap = "occ_map_buggy__.ppm";
		
		coop_ships_IDs.push_back("1");
		
		nbrAShips = 0;
		nbrAShipsMax = 10;
		nbrExtraObstacles = 10;
		
		distSH.push_back(2);
// 		distSH.push_back(2);

		distSO.push_back(2);
// 		distSO.push_back(2);

		distSA.push_back({2});
// 		distSA.push_back({2});

		distSS.push_back({2});

		speedbounds.push_back(2.0);
		
		Ninterp = 2;
		
		ratioToFirst.push_back(1.0);
// 		ratioToFirst.push_back(0.6);
		
		withCollReg = 0;
		ENDTIME = OPEN_WITH_FIXED_RATIOS;
		
		tEnd_nu.push_back(40);
// 		tEnd_nu.push_back(40);

		collavoidSH = 1;
		collavoidSO = 0;
		collavoidSA = 0;
		collavoidSS = 0;
		speedLimit  = 0;
		
		tWgt.push_back(100);
// 		tWgt.push_back(1000);

		eWgt.push_back(0.05);
// 		eWgt.push_back(0.5);

		fWgt.push_back({10000, 10000, nan(""), 500, nan(""), nan(""), nan(""), nan("")});
// 		fWgt.push_back({10000, 10000, nan(""), 500, nan(""), nan(""), nan(""), nan("")});
		
		cWgt.push_back({100, 0.2, 0.5});
// 		cWgt.push_back({100, 0.1, 0.5});

		plow_nu.push_back(10);
// 		plow_nu.push_back(10);
		pupp_nu.push_back(1000);
// 		pupp_nu.push_back(1000);
		
		scale_t.push_back(10);
// 		scale_t.push_back(100);
		
		scale.push_back({ 100., 100., 1, 10, 1, 1, 1, 100, 1. });
// 		scale.push_back({ 100., 100., 1, 10, 1, 1, 1, 100, 1. });
		
		scale_ctrl.push_back({ 0.01, 2 });
// 		scale_ctrl.push_back({ 0.01, 2 });

		viewer_t_scale = 10.0;
		
		stepSize = 15.0;
		
		buffer = 2.;
		ned_map_delta .push_back(2);
		ned_map_delta .push_back(2);
		ned_grid_delta.push_back(1);
		ned_grid_delta.push_back(1);
		
		soft_final_constr = 0;
	}
	
	string matlab_file = "temp1.m";
	
	int nbrCShips = coop_ships_IDs.size();

	vector <BaseModel*> ship_model ( nbrCShips );
	 
	for(int k=0; k<nbrCShips; k++) {
		if(vehicle == SHIP) {
			ship_model[k] = get_ship_model (coop_ships_IDs[k]);
		} else {
			ship_model[k] = get_buggy_model(coop_ships_IDs[k]);
		}
	}
	
	GeometricConstr geom_constr(nbrCShips, nbrAShips);
	geom_constr.lat0 = lat0;
	geom_constr.lon0 = lon0;
	geom_constr.hgt0 = hgt0;
	geom_constr.loadHarbPolygon(harb_f[0], grid_f[0]);
	geom_constr.loadGuidPolygon(guid_f[0]);
	for(int i=0; i<obst_f.size(); i++) geom_constr.loadObstPolygon(obst_f[i], obst_grid_f[i]);
	for(int i=0; i<guid_f.size(); i++) geom_constr.loadGuidPolygon(guid_f[i]);
	for(int i=0; i<spee_f.size(); i++) geom_constr.loadSpeePolygon(spee_f[i]);

	if(vehicle == SHIP) {
		for ( int i=0; i<nbrCShips; i++  ) geom_constr.initializeCShip   (i, ship_model[i]->L, ship_model[i]->BOA);
	} else {
		for ( int i=0; i<nbrCShips; i++  ) geom_constr.initializeCShip   (i, ship_model[i]->L, 0.75*ship_model[i]->L);
	}

	nbrObstacles = geom_constr.getNbrObstacles();
	nbrObstacles_init = nbrObstacles;
	nbrObstaclesMax = nbrObstacles_init + nbrExtraObstacles;
	geom_constr.transformHarbPolygonsToNED();
	geom_constr.transformObstPolygonsToNED();
	geom_constr.transformGuidPolygonsToNED();
	geom_constr.transformSpeePolygonsToNED();

	geom_constr.fOccMap = fOccMap;
	geom_constr.tOccMap = tOccMap;

	vector<vector<Vec2d>> simple_path_solutions(nbrCShips);

	string viewerID = "Optimal Trajectories";
	bool showDF = 1;
	bool showDG = 1;
	bool showHM = 1;
	
	vector<double> tEnd(tEnd_nu);
	vector<double> plow(plow_nu);
	vector<double> pupp(pupp_nu);
	for(int k=0; k<nbrCShips; k++) {
		if(ENDTIME == COMMON_AND_FIXED      ) tEnd[k] = tEnd_nu[k]/scale_t[0];
		if(ENDTIME == OPEN_WITH_FIXED_RATIOS) tEnd[k] = tEnd_nu[k]/scale_t[0];
		if(ENDTIME == OPEN                  ) tEnd[k] = tEnd_nu[k]/scale_t[k];
		plow[k] = plow_nu[k]/scale_t[k];
		pupp[k] = pupp_nu[k]/scale_t[k];
	}
	
	vector<int> NmodelS(nbrCShips);
	vector<int> NmodelC(nbrCShips);
	vector<int> NS(nbrCShips+1);
	vector<int> NC(nbrCShips+1);
	NS[0] = 0;
	NC[0] = 0;

	for ( int i=0; i<nbrCShips; i++ ) {
		NmodelS[i] = ship_model[i]->nbrS;
		NmodelC[i] = ship_model[i]->nbrC;
		NS[i+1] = NS[i] + NmodelS[i] + NmodelC[i] + 1;
		NC[i+1] = NC[i] + NmodelC[i];
	}
	
	
	tw::TWparameter twparameter ( transxml );
	map<string,string> args = twparameter.Arguments ( argv,argc );
	XMLNode *xml_shipviewer = tw::TWparameter::ReadParams ( shipvxml );
	
	tw::TWdimension TWdim;
	TWdim.ID     = viewerID;
	TWdim.n_dis  = twparameter.NDIS;
	TWdim.n_ode  = NS[nbrCShips];
	TWdim.n_ctrl = NC[nbrCShips];
	
	if(ENDTIME == COMMON_AND_FIXED      ) TWdim.n_param = 0;
	if(ENDTIME == OPEN_WITH_FIXED_RATIOS) TWdim.n_param = 1;
	if(ENDTIME == OPEN                  ) TWdim.n_param = nbrCShips;
	
	TWdim.n_rand  =              collavoidSH   *nbrCShips             *((TWdim.n_dis-1)*(Ninterp+1)+1)
					+            collavoidSO   *nbrCShips*nbrObstaclesMax*((TWdim.n_dis-1)*(Ninterp+1)+1)
					+            collavoidSA   *nbrCShips*nbrAShipsMax   *((TWdim.n_dis-1)*(Ninterp+1)+1);
	if(ENDTIME == COMMON_AND_FIXED      ) TWdim.n_rand  += collavoidSS*nChoosek(nbrCShips,2)*((TWdim.n_dis-1)*(Ninterp+1)+1);
	if(ENDTIME == OPEN_WITH_FIXED_RATIOS) TWdim.n_rand  += collavoidSS*nChoosek(nbrCShips,2)*((TWdim.n_dis-1)*(Ninterp+1)+1);
	TWdim.n_rand  += speedLimit*nbrCShips*TWdim.n_dis;
	
	
	twparameter.showDF = showDF;
	twparameter.showDG = showDG;
	twparameter.showHM = showHM; 
	
	Viewer *viewer = 0;
	if ( twparameter.PLOT ) viewer = new tw::Viewer ( &twparameter );
	
	tw::TWfolder folder ( &twparameter,0 );
	
	OptTraj ph ( xml_shipviewer, TWdim, nbrCShips, nbrAShips, nbrAShipsMax, nbrObstacles, nbrObstacles_init, nbrObstaclesMax, ship_model, geom_constr, 
				distSH, distSO, distSA, distSS, speedbounds,
				ENDTIME, tEnd, 
				collavoidSH, collavoidSO, collavoidSA, collavoidSS, speedLimit,
				Ninterp,  
				NS, NC, NmodelS, NmodelC, scale, scale_ctrl, scale_t, viewer_t_scale,
				tWgt, eWgt, fWgt, cWgt, plow, pupp, 
				ratioToFirst, withCollReg, soft_final_constr, &folder );
	ph0 = &ph;
	folder0 = &folder;
	
	ph.setSolver ( &twparameter );
	folder.Add ( &ph );
	folder.Init();
	
	uint8_t packetID;

	udp_client_server::udp_server server(ipad, port2);

	int idx_viewer = 0;
	idx_viewer0 = &idx_viewer;

	int solSize;
	vector<int> cur_idx(nbrCShips);
	for(int k=0; k<nbrCShips; k++) cur_idx[k] = 0;
	
	bool cs_first_time = true;
	bool cs_init = false;
	bool as_init = false;
	bool os_init = false;
	uint32_t cs_time = 1;
	uint32_t as_time = 2;
	uint32_t os_time = 3;
	
	vector<bool> fs_init(nbrCShips);
	for(int k=0; k<nbrCShips; k++) fs_init[k] = false;
	vector<bool> fs_updated(nbrCShips);
	for(int k=0; k<nbrCShips; k++) fs_updated[k] = false;
	bool first_time = true;
	
	
	double x0, y0, z0;
	while(1) {

		char recv_msg[MaxsizeofMsg];
		int sizeofMsg = server.timed_recv(recv_msg, MaxsizeofMsg, 2000);

		if(sizeofMsg == -1) {
			log("timed out or error occured in \"int udp_server::timed_recv(char *msg, size_t max_size)\"");
			cout << "sizeofMsg: " << sizeofMsg << endl;
		} else {
			memcpy(&packetID, recv_msg+2, 1); 

			if(packetID == DP_0XA0) {
				cout << "DP_0XA0" << endl;
				dp_0xA0 datadict;
				datadict.from_msg(recv_msg);
			}
			if(packetID == DP_0XA1) {
				cout << "DP_0XA1" << endl;
				dp_0xA1 datadict;
				datadict.from_msg(recv_msg);
			}
			if(packetID == DP_0XA4) {
				cout << "DP_0XA4" << endl;
				dp_0xA4 datadict;
				datadict.from_msg(recv_msg);
			}
			if(packetID == DP_0XAC) {
				cout << "DP_0XAC" << endl;
				dp_0xAC datadict;
				datadict.from_msg(recv_msg);
				datadict.toTW(lat0, lon0, hgt0);
				cs_init = true;
				cs_time = datadict.header[3].value;
			}
			if(packetID == DP_0XAD) {
				cout << "DP_0XAD"   << endl;
				dp_0xAD datadict;
				datadict.from_msg(recv_msg);
				datadict.toTW(lat0, lon0, hgt0);
				as_init = true;
				as_time = datadict.header[3].value;
			}
			if(packetID == DP_0XAE) {
				cout << "DP_0XAE"   << endl;
				dp_0xAE datadict;
				datadict.from_msg(recv_msg);
				datadict.toTW();
				os_init = true;
				os_time = datadict.header[3].value;
			}
			if(packetID == DP_0XAF) {
				cout << "DP_0XAF" << endl;
				uint16_t prodID;
				memcpy(&prodID, recv_msg+3, 2);
				dp_0xAF datadict;
				datadict.from_msg(recv_msg);
				datadict.toTW(prodIDmap[prodID], lat0, lon0, hgt0);
				fs_init[prodIDmap[prodID]] = true;
				fs_updated[prodIDmap[prodID]] = true;
			}

			if((cs_time==as_time) && (as_time==os_time) && cs_init && as_init && os_init && is_all_true(fs_init)) {
				cs_init = false;
				as_init = false;
				os_init = false;
				
				if(cs_first_time) {
					for(int k=0; k<nbrCShips; k++) {
						ph.var_b_final[k][0] = ph.var_b_init[k][0];
						ph.var_b_final[k][1] = ph.var_b_init[k][1];
					}
				}
				cs_first_time = false;
				
				if(!is_all_false(fs_updated)) {
					// Simple Path Search
					Eigen::Vector2d d_xy_ned_map(ned_map_delta[0], ned_map_delta[1]);
					Eigen::Vector2d d_xy_ned_grid(ned_grid_delta[0], ned_grid_delta[1]);
					geom_constr.initializeOccMap(d_xy_ned_map, buffer);
					
					for(int k=0; k<nbrCShips; k++) {
						if(fs_updated[k]) {
							double lat, lon, hgt;
							NEDtoLLA(ph.var_b_final[k][0]*ph.scale[k][0], ph.var_b_final[k][1]*ph.scale[k][1], 0, lat0, lon0, hgt0, &lat, &lon, &hgt);
							Eigen::Vector2d start_lla(lat, lon);
							NEDtoLLA(ph.final_s[k][0]*ph.scale[k][0], ph.final_s[k][1]*ph.scale[k][1], 0, lat0, lon0, hgt0, &lat, &lon, &hgt);
							Eigen::Vector2d goal_lla (lat, lon);
							dsl::SimplePathSearch simple_path_search(geom_constr, d_xy_ned_grid, start_lla, goal_lla);
							simple_path_solutions[k] = simple_path_search.solution;
							
							cur_idx[k] = 0;
							fs_updated[k] = false;
						}
					}
					// End of Simple Path Search
				}
				
				//update var_b_final
				for(int k=0; k<nbrCShips; k++) {
					solSize = simple_path_solutions[k].size();
					
					LLAtoNED(simple_path_solutions[k][cur_idx[k]].x, 
									simple_path_solutions[k][cur_idx[k]].y, 0, lat0, lon0, hgt0, &x0, &y0, &z0);

					while(  (ph.var_b_final[k][0]*scale[k][0]-x0)*(ph.var_b_final[k][0]*scale[k][0]-x0) + 
							(ph.var_b_final[k][1]*scale[k][1]-y0)*(ph.var_b_final[k][1]*scale[k][1]-y0) < stepSize*stepSize &&
							cur_idx[k] < solSize-1 ) {
						cur_idx[k]++;
						LLAtoNED(simple_path_solutions[k][cur_idx[k]].x, 
									simple_path_solutions[k][cur_idx[k]].y, 0, lat0, lon0, hgt0, &x0, &y0, &z0);
					}
					
					ph.var_b_final[k][0] = x0/scale[k][0];
					ph.var_b_final[k][1] = y0/scale[k][1];
// 					for(int i=2; i<9; i++) ph.var_b_final[k][i] = nan("");
// 					ph.var_b_final[k][3] = 0;
				}
				
				for(int k=0; k<ph.nbrCShips; k++) {
					if(ph.solver->X[ph.x_index ( 0,ph.NS[k]+2 )]*ph.scale[k][2]*180./M_PI - ph.var_b_init[k][2]*ph.scale[k][2]*180./M_PI > 180) {
						ph.var_b_init[k][2] = (ph.var_b_init[k][2]*ph.scale[k][2] + 2*M_PI)/ph.scale[k][2];
					}
					if(ph.solver->X[ph.x_index ( 0,ph.NS[k]+2 )]*ph.scale[k][2]*180./M_PI - ph.var_b_init[k][2]*ph.scale[k][2]*180./M_PI < -180) {
						ph.var_b_init[k][2] = (ph.var_b_init[k][2]*ph.scale[k][2] - 2*M_PI)/ph.scale[k][2];
					}
				}
				
				//Course Correction
// 				if(!first_time) {
// 					for(int k=0; k<nbrCShips; k++) {
// 						vector<double> pertubation(ph.n_ode);
// 						for (int i = 0; i < 6       ; i++) {
// 							pertubation[i] = ph.var_b_init[k][i] - ph.solver->X[ph.x_index ( 0,ph.NS[k]+i )];
// 							cout << "perturbation  " << pertubation[i] << endl;
// 						}
// 						for (int i = 6; i < ph.n_ode; i++) {
// 							pertubation[i] = 0.0;
// 						}
// 						for(int i=0; i<20; i++) {
// 							cout << ph.solver->X[ph.x_index ( i,ph.NS[k]+0 )]*ph.scale[k][0] << "    ";
// 							cout << ph.solver->X[ph.x_index ( i,ph.NS[k]+1 )]*ph.scale[k][1] << "    ";
// 							cout << ph.solver->X[ph.x_index ( i,ph.NS[k]+2 )]*ph.scale[k][2] << "    ";
// 							cout << endl;
// 						}
// 						korrektur_startposition(pertubation);
// 						for(int i=0; i<20; i++) {
// 							cout << ph.solver->X[ph.x_index ( i,ph.NS[k]+0 )]*ph.scale[k][0] << "    ";
// 							cout << ph.solver->X[ph.x_index ( i,ph.NS[k]+1 )]*ph.scale[k][1] << "    ";
// 							cout << ph.solver->X[ph.x_index ( i,ph.NS[k]+2 )]*ph.scale[k][2] << "    ";
// 							cout << endl;
// 						}
// 						cout << endl;
// 						ph.ToMATLAB ( matlab_file );
// 					}
// 				}
				
				
				folder.Init(); //This sets X to zero
				
				for(int k=0; k<nbrCShips; k++) {
					for ( int i=0; i<NS[k+1]-NS[k]; i++ ) ph.solver->X[ph.x_index ( 0,NS[k]+i )] = ph.var_b_init[k][i];
				}

				if(first_time) {
					for(int k=0; k<nbrCShips; k++) {
						for ( int i=0; i<NC[k+1]-NC[k]; i++ ) {
							for(int j=0; j<ph.n_dis; j++) ph.solver->X[ph.u_index ( j,NC[k]+i )] = 0;
						}
						ph.solver->X[ph.x_index ( 0,NS[k]+6 )] = 0;
						if(vehicle == SHIP ) ph.solver->X[ph.x_index ( 0,NS[k]+7 )] = 0.1/ph.scale[k][7];
						if(vehicle == BUGGY) ph.solver->X[ph.x_index ( 0,NS[k]+7 )] = 19/ph.scale[k][7];
					}
					ph.solver->Integrate ( twparameter.butchertableau );
					first_time = false;
				} else {
					ph.FromMATLAB ( matlab_file );
				}
				
				
				
				if(idx_viewer > 15) folder.Init ( viewer );
				folder.Loop(0);
				if (viewer) SetFloatTime(viewer, viewer_t_scale);
				if(idx_viewer > 14) folder.Init ( viewer );
				idx_viewer++;
				folder.Loop(0);
				ph.ToMATLAB ( matlab_file );
				
				
				for(int k=0; k<nbrCShips; k++) {
					dp_base* datadict;
					if(vehicle == SHIP ) datadict = new dp_0xA3;
					if(vehicle == BUGGY) datadict = new dp_0xA2;
					datadict->header[0].value = 1;
					if(vehicle == SHIP ) datadict->header[1].value = DP_0XA3;
					if(vehicle == BUGGY) datadict->header[1].value = DP_0XA2;
					datadict->header[2].value = k;
					datadict->header[3].value = 0;
					datadict->header[4].value = 0;
					datadict->header[5].value = 0;
					datadict->footer[0].value = 0;
					if(vehicle == SHIP ) datadict->fromTW(k, lat0, lon0, hgt0, cs_time, 5);
					if(vehicle == BUGGY) datadict->fromTW(k, lat0, lon0, hgt0, cs_time, 1);
					std::this_thread::sleep_for (std::chrono::milliseconds(100));
					cout << "################## Sending message" << "  " << time(NULL) << endl;
					datadict->udp_send(ipad, port1);
					cout << "################## Message sent" << "  " << time(NULL) << endl;
					delete datadict;
				}
				
// 				for(int k=0; k<nbrCShips; k++) {
// 					dp_base* datadict;
// 					if(vehicle == SHIP ) datadict = new dp_0xA3;
// 					if(vehicle == BUGGY) datadict = new dp_0xA2;
// 					datadict->header[0].value = 21845;
// 					if(vehicle == SHIP ) datadict->header[1].value = DP_0XA3;
// 					if(vehicle == BUGGY) datadict->header[1].value = DP_0XA2;
// 					datadict->header[2].value = k+5;
// 					datadict->header[3].value = 0;
// 					datadict->header[4].value = 0;
// 					datadict->header[5].value = 0;
// 					datadict->footer[0].value = 0;
// 					if(vehicle == SHIP ) datadict->fromTW(k, lat0, lon0, hgt0, cs_time, 5);
// 					if(vehicle == BUGGY) datadict->fromTW(k, lat0, lon0, hgt0, cs_time, 0.1);
// 					cout << *datadict << endl;
// 					exit(0);
// 					std::this_thread::sleep_for (std::chrono::milliseconds(100));
// 					cout << "################## Sending message" << "  " << time(NULL) << endl;
// 					datadict->udp_send("10.121.1.2", 20024);
// 					cout << "################## Message sent" << "  " << time(NULL) << endl;
// 					delete datadict;
// 				}
			} 
		}
	}

	delete viewer;
	delete xml_shipviewer;
	for ( int i=0; i<nbrCShips; i++ ) delete ship_model[i];

	return 0;
}



