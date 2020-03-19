/*    Optimal Trajectory Calculations     */

#ifdef WIN32
#include "windows.h"
#endif
#include "optimal_traj_Nships.h"
#include "gridpathsearch.h"
#include "global.h"

#include <fstream>
#include <time.h>

 
void log( const std::string &text )
{
    std::ofstream log_file(
        "log_file.txt", std::ios_base::out | std::ios_base::app );
    log_file << text << endl;
}


int main ( int argv, char* argc[] ) {
	
	string transxml = "transworhp.xml";
	string shipvxml = "shipviewer.xml";
	vector<string> harb_f   = {"RostockHarbor/harbor_simplified.txt"};
	vector<string> grid_f   = {"RostockHarbor/harbor_simplified.xml"};
	vector<string> guid_f   = {"guide.txt"};
	vector<string> obst_f   = {"RostockHarbor/obstacle_01_simplified.txt",
							   "RostockHarbor/obstacle_02_simplified.txt",
							   "RostockHarbor/obstacle_03_simplified.txt",
							   "RostockHarbor/obstacle_04_simplified.txt",
							   "RostockHarbor/obstacle_05_simplified.txt",
							   "RostockHarbor/obstacle_06_simplified.txt",
							   "RostockHarbor/obstacle_07_simplified.txt",
							   "RostockHarbor/obstacle_08_simplified.txt",
							   "RostockHarbor/obstacle_09_simplified.txt",
							   "RostockHarbor/obstacle_10_simplified.txt"
	};
	vector<string> obst_grid_f   = {""/*"RostockHarbor/obstacle_01_simplified.xml"*/,
								    ""/*"RostockHarbor/obstacle_02_simplified.xml"*/,
									""/*"RostockHarbor/obstacle_03_simplified.xml"*/,
									""/*"RostockHarbor/obstacle_04_simplified.xml"*/,
									""/*"RostockHarbor/obstacle_05_simplified.xml"*/,
									""/*"RostockHarbor/obstacle_06_simplified.xml"*/,
									""/*"RostockHarbor/obstacle_07_simplified.xml"*/,
									""/*"RostockHarbor/obstacle_08_simplified.xml"*/,
									""/*"RostockHarbor/obstacle_09_simplified.xml"*/,
									""/*"RostockHarbor/obstacle_10_simplified.xml"*/,
	};
	vector<string> spee_f   = {"speedbound.txt"};
	
	string ipad = "127.0.0.1";
	int    port = 1502;
	
	double lat0 = 54.17057475;
	double lon0 = 12.10074142;
	double hgt0 = 0;
	
	string fMatlab;
	string tMatlab  = "temp1.m";
	
	vector <string> coop_ships_IDs = {"1","2"};
	
	vector <vector<double>> AIS_ships_dims = {
		{200., 47.5},
	};
	
	int nbrCShips = coop_ships_IDs.size();
	int nbrAShips = AIS_ships_dims.size();
	int nbrObstacles;
	
	vector <BaseModel*> ship_model ( nbrCShips );
	for(int i=0; i<nbrCShips; i++) ship_model[i] = get_ship_model(coop_ships_IDs[i]);
	
	GeometricConstr geom_constr(nbrCShips, nbrAShips);
	geom_constr.lat0 = lat0;
	geom_constr.lon0 = lon0;
	geom_constr.hgt0 = hgt0;
	geom_constr.loadHarbPolygon(harb_f[0], grid_f[0]);
	geom_constr.loadGuidPolygon(guid_f[0]);
	for(int i=0; i<obst_f.size(); i++) geom_constr.loadObstPolygon(obst_f[i], obst_grid_f[i]);
	geom_constr.loadSpeePolygon(spee_f[0]);

	for ( int i=0; i<nbrCShips; i++  ) geom_constr.initializeCShip   (i, ship_model[i]->L, ship_model[i]->BOA);
	for ( int i=0; i<nbrAShips; i++  ) geom_constr.initializeAShip(i, AIS_ships_dims[i][0], AIS_ships_dims[i][1]);
	
	nbrObstacles = geom_constr.getNbrObstacles();
	geom_constr.transformHarbPolygonsToNED();
	geom_constr.transformObstPolygonsToNED();
	geom_constr.transformGuidPolygonsToNED();
	geom_constr.transformSpeePolygonsToNED();
	
	vector <vector<double>> initial_nu = {
		{    54.183615, 12.091236,  +160.0,   0,  0 ,   0 ,   0 ,  0. ,  0 },
		{    54.165546, 12.129062,  +270.0,   0,  0 ,   0 ,   0 ,  0. ,  0 },
	};
	// the final heading has to be specified for the arc-search
// 	vector <vector<double>> final_nu = {
// 		{    54.160538, 12.118279,   nan("") ,  0 , 0 ,  0 ,  nan("") ,  nan("") ,  nan("") },
// 		{    54.148667, 12.093325,   +2000  ,  0,  0 ,  0 ,  nan("") ,  nan("") ,  nan("") },.
// 	};
	vector <vector<double>> final_nu = {
		{    54.160538, 12.118279,   nan("") ,  0 , 0 ,  0 ,  nan("") ,  nan("") ,  nan("") },
		{    54.148667, 12.093325,   nan("") ,  0,  0 ,  0 ,  nan("") ,  nan("") ,  nan("") },
	};
	// vector <vector<double>> final_nu = {
	// 	{    54.160538, 12.118279,   nan("") ,  0 , 0 ,  0 ,  nan("") ,  nan("") ,  nan("") },
	// 	{    54.148667, 12.093325,   nan("") ,  0,  0 ,  0 ,  nan("") ,  nan("") ,  nan("") },
	// };
	
	vector<vector<Vec2d>> simple_path_solutions(nbrCShips);
	
	// Simple Path Search
	// Eigen::Vector2d d_xy_ned_map(20., 20.);
	// double buffer = 0.;
	// Eigen::Vector2d d_xy_ned_grid(10., 10.);
	// geom_constr.create_new_occ_map=false;
	// geom_constr.initializeOccMap(d_xy_ned_map, buffer);
	
	// for(int k=1; k<nbrCShips; k++) {
	// 	Eigen::Vector2d start_lla(initial_nu[k][0], initial_nu[k][1]);
	// 	Eigen::Vector2d goal_lla (final_nu  [k][0], final_nu  [k][1]);
	// 	dsl::SimplePathSearch simple_path_search(geom_constr, d_xy_ned_grid, start_lla, goal_lla);

	// 	simple_path_solutions[k] = simple_path_search.solution;
	// }
	// End of Simple Path Search

	// Arc Path Search

	vector<vector<Eigen::Vector3d>> arc_path_solutions(nbrCShips);

	Eigen::Vector2d d_xy_ned_map(20., 20.);
	double buffer = 0.;
	Eigen::Vector2d d_xy_ned_grid(10., 10.);
	geom_constr.create_new_occ_map=false;
	geom_constr.initializeOccMap(d_xy_ned_map, buffer);
	
	for(int k=0; k<nbrCShips; k++) {
		Eigen::Vector2d start_lla(initial_nu[k][0], initial_nu[k][1]);
		Eigen::Vector2d goal_lla (final_nu  [k][0], final_nu  [k][1]);
		dsl::SimplePathSearch simple_path_search(geom_constr, d_xy_ned_grid, start_lla, goal_lla);

		simple_path_solutions[k] = simple_path_search.solution;
	}
	// End of Simple Path Search

	// Arc Path Search

// 	vector<vector<Eigen::Vector3d>> arc_path_solutions(nbrCShips);
// 
// 	Eigen::Vector2d d_xy_ned_map(20., 20.);
// 	double buffer = 0.;
// 	Eigen::Vector3d d_nepsi_grid(10., 10., 2*M_PI/100.);
// 	geom_constr.create_new_occ_map=false;
// 	geom_constr.initializeOccMap(d_xy_ned_map, buffer);
// 	
// 	for(int k=1; k<nbrCShips; k++) {
// 		Eigen::Vector3d start_lla(initial_nu[k][0], initial_nu[k][1], initial_nu[k][2]/180.*M_PI );
// 		Eigen::Vector3d goal_lla (final_nu  [k][0], final_nu  [k][1], final_nu  [k][2]/180.*M_PI );
// 		dsl::ArcPathSearch arc_path_search(geom_constr, d_nepsi_grid, start_lla, goal_lla);
// 
// 		arc_path_solutions[k] = arc_path_search.solution;
// 	}
	// End of Arc Path Search

	exit(1);
	
	vector<double>         distSH = {50, 50};
	vector<double>         distSO = {60, 60};
	vector<vector<double>> distSA = {
		{100, 100}
	};
	vector<vector<double>> distSS = {
		{100}
	};
	vector<double> speedbounds = {2.0};
	
	int Ninterp = 2;
	
	vector <double> ratioToFirst = {1.0, 0.6};
	
	bool withCollReg = 0;
	EndTimeBool ENDTIME = OPEN_WITH_FIXED_RATIOS;
	vector<double> tEnd_nu = {300,300};
	bool collavoidSH = 1;
	bool collavoidSO = 1;
	bool collavoidSA = 0;
	bool collavoidSS = 1;
	bool speedLimit  = 1;
	
	vector<double>         tWgt = {2,2};
	vector<double>         eWgt = {1,1};
	vector<vector<double>> fWgt = {
								{200, 200, nan(""), 10, 10, 10, nan(""), nan("")},
								{200, 200, nan(""), 10, 10, 10, nan(""), nan("")},
								};
	vector<vector<double>> cWgt = {
								{10, 10, 0.003},
								{10, 10, 0.003},
								};
	vector<double>         plow_nu = {100,100};
	vector<double>         pupp_nu = {4000,4000};
	
	vector<double>          scale_t = {1500, 700};
	
	vector <vector<double>> scale = {
		{ 2000., 2000., 6.0, 3.5, 0.1, 0.001, 0.3, 1., 1. },
		{ 2000., 2000., 6.0, 3.5, 0.1, 0.001, 0.3, 1., 1. },
	};
	
	vector <vector<double>> scale_ctrl = {
		{ 0.001, 0.1 },
		{ 0.02, 0.1 },
	};
	
	
	double viewer_t_scale = 120.0;
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

	for ( int i=0; i<nbrCShips; i++  ) geom_constr.initializeCShip   (i, ship_model[i]->L, ship_model[i]->BOA);
	for ( int i=0; i<nbrAShips; i++  ) geom_constr.initializeAShip(i, AIS_ships_dims[i][0], AIS_ships_dims[i][1]);
	
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

	vector<vector<double>> initial(initial_nu);
	vector <vector<double>> sub_final_nu(final_nu);
	vector<vector<double>> final(sub_final_nu);
	
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
					+            collavoidSO   *nbrCShips*nbrObstacles*((TWdim.n_dis-1)*(Ninterp+1)+1)
					+            collavoidSA*nbrCShips*nbrAShips *((TWdim.n_dis-1)*(Ninterp+1)+1);
	if(ENDTIME == COMMON_AND_FIXED      ) TWdim.n_rand  += collavoidSS*nChoosek(nbrCShips,2)*((TWdim.n_dis-1)*(Ninterp+1)+1);
	if(ENDTIME == OPEN_WITH_FIXED_RATIOS) TWdim.n_rand  += collavoidSS*nChoosek(nbrCShips,2)*((TWdim.n_dis-1)*(Ninterp+1)+1);
	TWdim.n_rand  += speedLimit*nbrCShips*TWdim.n_dis;
	
	
	twparameter.showDF = showDF;
	twparameter.showDG = showDG;
	twparameter.showHM = showHM; 
	
	Viewer *viewer = 0;
	if ( twparameter.PLOT ) viewer = new tw::Viewer ( &twparameter );
	
	tw::TWfolder folder ( &twparameter,0 );
	
	int nbrAShipsMax = 10;
	int nbrObstacles_init  = nbrObstacles;
	int nbrObstaclesMax = nbrObstacles;
	
	OptTraj ph ( xml_shipviewer, TWdim, nbrCShips, nbrAShips, nbrAShipsMax, nbrObstacles, nbrObstacles_init, nbrObstaclesMax, ship_model, geom_constr, 
				distSH, distSO, distSA, distSS, speedbounds,
				ENDTIME, tEnd, 
				collavoidSH, collavoidSO, collavoidSA, collavoidSS, speedLimit,
				Ninterp,  
				NS, NC, NmodelS, NmodelC, scale, scale_ctrl, scale_t, viewer_t_scale,
				tWgt, eWgt, fWgt, cWgt, plow, pupp, 
				ratioToFirst, withCollReg, &folder );
	ph0 = &ph;
	
	ph.setSolver ( &twparameter );
	folder.Add ( &ph );
	folder.Init();
	
	int nbrSubCalc = 10;
	
	for(int index=0; index<nbrSubCalc-1; index++) { 
		
		if(index==0) fMatlab  = "";
		if(index!=0) fMatlab  = "temp1.m";
		
		int solSize;
		for(int k=0; k<nbrCShips; k++) {
			solSize = simple_path_solutions[k].size();
			sub_final_nu[k][0] = simple_path_solutions[k][(solSize-1)*(index+1)/(nbrSubCalc-1)].x;
			sub_final_nu[k][1] = simple_path_solutions[k][(solSize-1)*(index+1)/(nbrSubCalc-1)].y;
		}

		double z0;
		for(int k=0; k<nbrCShips; k++) {
			LLAtoNED(initial_nu[k][0], initial_nu[k][1], 0, lat0, lon0, hgt0, &(initial[k][0]), &(initial[k][1]), &z0);
			LLAtoNED(sub_final_nu  [k][0], sub_final_nu  [k][1], 0, lat0, lon0, hgt0, &(final  [k][0]), &(final  [k][1]), &z0);
			initial[k][2] = initial_nu[k][2]*M_PI/180;
			final  [k][2] = sub_final_nu  [k][2]*M_PI/180;
		}
		vector<vector<double>> initial_guess(initial);
		for(int k=0; k<nbrCShips; k++) initial_guess[k][7] = 1;
		
		for(int k=0; k<nbrCShips; k++) {
			for(int i=0; i<NmodelS[k]+NmodelC[k]+1; i++) {
				ph.var_b_init [k][i] = initial[k][i]/scale[k][i];
				ph.var_b_final[k][i] = final[k][i]/scale[k][i];
				initial_guess[k][i] /= scale[k][i];
			}
		}

		
			folder.Init();

			if(fMatlab.empty()) {
				for ( int k=0; k<nbrCShips; k++ ) {
					for ( int i=0; i<NmodelS[k]+NmodelC[k]+1; i++ ) ph.solver->X[ph.x_index ( 0,NS[k]+i )] = initial_guess[k][i];
				}
				int steps = ph.solver->Integrate ( twparameter.butchertableau );
				cout << "Integrationsschritte: " << steps << endl;
			}
			if(!fMatlab.empty()) ph.FromMATLAB ( fMatlab );
			
// 			if (viewer) SetFloatTime(viewer, 0.1*viewer_t_scale);
// 			folder.Init ( viewer );
			folder.Loop(0);
			if (viewer) SetFloatTime(viewer, viewer_t_scale);
			/*if(index==nbrSubCalc-2)*/ folder.Init ( viewer );
			folder.Loop(0);
			if(!tMatlab.empty()) ph.ToMATLAB ( tMatlab );

	}
	delete viewer;
	delete xml_shipviewer;
	for ( int i=0; i<nbrCShips; i++ ) delete ship_model[i];

	return 0;
}



