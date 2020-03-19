/*    Optimal Trajectory Calculations     */

#ifdef WIN32
#include "windows.h"
#endif
#include "optimal_traj_Nships.h"
#include "global.h"

#include <fstream>
#include <time.h>

tw::TWfolder *myfolder;
int main ( int argv, char* argc[] ) {

	
	if (argv != 10){
		printf("Number of input arguments: %d\n",argv);
		printf("Please follow the following input format:\n");
		printf("starting position: (lan lon psi) end position: (lat lon)\n");
		return 1;
	}

    double start_lat = atof(argc[4]);
    double start_lon = atof(argc[5]);
	double start_psi = atof(argc[6]);
	double end_lat = atof(argc[7]);
    double end_lon = atof(argc[8]);
	string resultfile = argc[9];

	string transxml = "transworhp.xml";
	string shipvxml = "shipviewer_motion_planning.xml";
	
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
	
	double lat0 = 54.17057475;
	double lon0 = 12.10074142;
	double hgt0 = 0;
	
	//vector<string> harb_f = {"RostockHarbor/harbor_simplified.txt"};
	//vector<string> grid_f;
	//vector<string> guid_f;
	//vector<string> obst_f;
	//vector<string> obst_grid_f;
	//vector<string> spee_f;
	
	//double lat0;
	//double lon0;
	//double hgt0;
	
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
	bool collavoidSH = 1;
	bool collavoidSO = 0;
	bool collavoidSA = 0;
	bool collavoidSS = 0;
	bool speedLimit = 0;
	
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
	
	bool soft_final_constr;
	
	{
	
		
		coop_ships_IDs.push_back("1");
		
		nbrAShips = 0;
		nbrAShipsMax = 10;
		nbrExtraObstacles = 10;
		
		distSH.push_back(10);
		distSO.push_back(10);
		distSA.push_back({100});
		distSS.push_back({100});
		speedbounds.push_back(2.0);

		
		ratioToFirst.push_back(1.0);
		
		withCollReg = 0;
		ENDTIME = OPEN_WITH_FIXED_RATIOS;
		tEnd_nu.push_back(200);

		
		//tWgt.push_back(200);
		//eWgt.push_back(2);
		tWgt.push_back(500);
		eWgt.push_back(0.2);
		fWgt.push_back({200, 200, nan(""), 10, 10, 10, nan(""), nan("")});
		cWgt.push_back({10, 10, 0.003});
		plow_nu.push_back(100);
		pupp_nu.push_back(5000);
		scale_t.push_back(1500);
		scale.push_back({ 2000., 2000., 6.0, 3.5, 0.1, 0.001, 0.3, 1., 1. });
		scale_ctrl.push_back({ 0.001, 0.1 });
		
		viewer_t_scale = 120.0;
		
		soft_final_constr = 0;
		
		int nbrCShips = coop_ships_IDs.size();
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
		//for ( int i=0; i<nbrAShips; i++  ) geom_constr.initializeAShip(i, AIS_ships_dims[i][0], AIS_ships_dims[i][1]);
	
		nbrObstacles = geom_constr.getNbrObstacles();
		geom_constr.transformHarbPolygonsToNED();
		geom_constr.transformObstPolygonsToNED();
		geom_constr.transformGuidPolygonsToNED();
		geom_constr.transformSpeePolygonsToNED();
		
		//GeometricConstr geom_constr(nbrCShips, nbrAShips);
		//geom_constr.loadHarbPolygon(harb_f[0], grid_f[0]);
		string viewerID = "Optimal Trajectories";
		
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
		
		vector <vector<double>> sub_final_nu(nbrCShips);
		for(int k=0; k<nbrCShips; k++) sub_final_nu[k].resize(3);

		
		tw::TWparameter twparameter ( transxml );
		map<string,string> args = twparameter.Arguments ( argv,&argc[0] );
		auto xml_shipviewer = tw::TWparameter::ReadParams ( shipvxml );
		
		tw::TWdimension TWdim;
		TWdim.ID     = viewerID;
		TWdim.n_dis  = twparameter.NDIS;
		TWdim.n_ode  = NS[nbrCShips];
		TWdim.n_ctrl = NC[nbrCShips];
		
		if(ENDTIME == COMMON_AND_FIXED      ) TWdim.n_param = 0;
		if(ENDTIME == OPEN_WITH_FIXED_RATIOS) TWdim.n_param = 1;
		if(ENDTIME == OPEN                  ) TWdim.n_param = nbrCShips;
		
		//TWdim.n_rand  = 0;
		TWdim.n_rand  =              collavoidSH   *nbrCShips             *((TWdim.n_dis-1)*(Ninterp+1)+1)
					+            collavoidSO   *nbrCShips*nbrObstacles*((TWdim.n_dis-1)*(Ninterp+1)+1)
					+            collavoidSA*nbrCShips*nbrAShips *((TWdim.n_dis-1)*(Ninterp+1)+1);

		if(ENDTIME == COMMON_AND_FIXED      ) TWdim.n_rand  += collavoidSS*nChoosek(nbrCShips,2)*((TWdim.n_dis-1)*(Ninterp+1)+1);
		if(ENDTIME == OPEN_WITH_FIXED_RATIOS) TWdim.n_rand  += collavoidSS*nChoosek(nbrCShips,2)*((TWdim.n_dis-1)*(Ninterp+1)+1);
		TWdim.n_rand  += speedLimit*nbrCShips*TWdim.n_dis;
		//TWdim.n_rand  = 0;

		tw::Viewer *viewer = 0;
		if ( twparameter.PLOT ) viewer = new tw::Viewer ( &twparameter );
		
		tw::TWfolder folder ( &twparameter,0 );
		OptTraj ph ( xml_shipviewer.get(), TWdim, nbrCShips, nbrAShips, nbrAShipsMax, nbrObstacles, nbrObstacles_init, nbrObstaclesMax, ship_model, geom_constr, 
					distSH, distSO, distSA, distSS, speedbounds,
					ENDTIME, tEnd, 
					collavoidSH, collavoidSO, collavoidSA, collavoidSS, speedLimit,
					Ninterp,  
					NS, NC, NmodelS, NmodelC, scale, scale_ctrl, scale_t, viewer_t_scale,
					tWgt, eWgt, fWgt, cWgt, plow, pupp, 
					ratioToFirst, withCollReg, soft_final_constr, &folder, resultfile);
		ph0 = &ph;
		folder0 = &folder;
		myfolder = & folder;

		ph.setSolver ( &twparameter );
		folder.Add ( &ph );
		folder.Init();

		double x0, y0, z0;
		double lat;
		double lon;
		double psi;
		double u;
		double v;
		double r;
		double delta;
		double n;
		
		vector <vector<double>> initial_nu;
		vector <vector<double>> final_nu;

		//create motion primitives
		//initial_nu.push_back({ 54.17057475, 12.10074142, 0.0, 3.0, 0, 0, 0, 36.14521, 0 });
		//final_nu  .push_back({  nan("") ,nan("") , -30.0  , 3.0 , 0 ,  0 , 0,  36.14521 ,  nan("") });

		//reoptimization
		//intial:
		//rechts vom eingang
		//initial_nu.push_back({ 54.184672, 12.111759, 340.0, 3.0, 0, 0, 0, 36.14521, 0 });
		//links vom eingang
		//initial_nu.push_back({ 54.198314, 12.072334, 170.0, 3.0, 0, 0, 0, 36.14521, 0 });
		
		initial_nu.push_back({ start_lat, start_lon, start_psi, 3.0, 0, 0, 0, 36.14521, 0 });
	
		//main that takes traj as input
		//initial_nu.push_back({ start_lat, start_lon, start_psi, 3.0, 0, 0, 0, 36.14521, 0 });
		
		//final:

		//main that takes traj as input
		//final_nu  .push_back({ end_lat,end_lon, nan(""), nan(""), nan(""), nan(""), nan(""), nan(""), nan("") });
		
		// halbe Strecke bis Rostock (MOTION_PRIMS_3)
		//final_nu  .push_back({ 54.144776, 12.092533, nan(""), nan(""), nan(""), nan(""), nan(""), nan(""), nan("") });
		
		//volle Strecke bis Rostock
		//final_nu  .push_back({ 54.094846, 12.121948, nan(""),  nan(""), nan(""), nan(""), nan(""), nan(""), nan("") });
		
		// bis ins park dock
		
		final_nu  .push_back({  end_lat, end_lon,   nan("") ,  nan(""), nan("") , nan("") ,  nan("") ,  nan("") ,  nan("") });

		lat   = initial_nu[0][0];
		lon   = initial_nu[0][1];
		psi   = initial_nu[0][2]*M_PI/180;
		u     = initial_nu[0][3];
		v     = initial_nu[0][4];
		r     = initial_nu[0][5]*M_PI/180;
		delta = initial_nu[0][6]*M_PI/180;
		n     = initial_nu[0][7]/60;
		
		
		LLAtoNED(lat, lon, 0, lat0, lon0, 0, &x0, &y0, &z0);
		ph0->var_b_init[0][0] = x0;
		ph0->var_b_init[0][1] = y0;
		ph0->var_b_init[0][2] = psi;
		ph0->var_b_init[0][3] = u;
		ph0->var_b_init[0][4] = v;
		ph0->var_b_init[0][5] = r;
		ph0->var_b_init[0][6] = delta;
		ph0->var_b_init[0][7] = n;
		for(int i=0; i<8; i++) ph0->var_b_init[0][i] /= ph0->scale[0][i];

		ph0->solver->X[ph0->x_index ( 0,ph0->NS[0]+0 )] = x0;
		ph0->solver->X[ph0->x_index ( 0,ph0->NS[0]+1 )] = y0;
		ph0->solver->X[ph0->x_index ( 0,ph0->NS[0]+2 )] = psi;
		ph0->solver->X[ph0->x_index ( 0,ph0->NS[0]+3 )] = u;
		ph0->solver->X[ph0->x_index ( 0,ph0->NS[0]+4 )] = v;
		ph0->solver->X[ph0->x_index ( 0,ph0->NS[0]+5 )] = r;
		ph0->solver->X[ph0->x_index ( 0,ph0->NS[0]+6 )] = delta;
		ph0->solver->X[ph0->x_index ( 0,ph0->NS[0]+7 )] = n;
		for(int i=0; i<8; i++) ph0->solver->X[ph0->x_index ( 0,ph0->NS[0]+i )] /= ph0->scale[0][i];
		
		
		
		lat   = final_nu[0][0];
		lon   = final_nu[0][1];
		psi   = final_nu[0][2]*M_PI/180;
		u     = final_nu[0][3];
		v     = final_nu[0][4];
		r     = final_nu[0][5]*M_PI/180;
		delta = final_nu[0][6]*M_PI/180;
		n     = final_nu[0][7]/60;
		
		LLAtoNED(lat, lon, 0, lat0, lon0, 0, &x0, &y0, &z0);
		ph0->var_b_final[0][0] = x0; 
		ph0->var_b_final[0][1] = y0; 
		ph0->var_b_final[0][2] = psi; 
		ph0->var_b_final[0][3] = u; 
		ph0->var_b_final[0][4] = v; 
		ph0->var_b_final[0][5] = r; 
		ph0->var_b_final[0][6] = delta; 
		ph0->var_b_final[0][7] = n; 
		ph0->var_b_final[0][8] = nan(""); 
		for(int i=0; i<9; i++) ph0->var_b_final[0][i] /= ph0->scale[0][i];

		
		
		folder.Init();

		string fMatlab  = "";
		if(!fMatlab.empty()) {
			ph.FromMATLAB ( fMatlab );
		} else {
			for ( int i=0; i<NC[1]-NC[0]; i++ ) {
				for(int j=0; j<ph.n_dis; j++) ph.solver->X[ph.u_index ( j,NC[0]+i )] = 0;
			}
			ph.solver->X[ph.x_index ( 0,NS[0]+6 )] = 0;
			ph.solver->X[ph.x_index ( 0,NS[0]+7 )] = 36.14521/60./ph.scale[0][7];
			ph.solver->Integrate ( twparameter.butchertableau );
		}
		

	 	if (viewer) SetFloatTime(viewer, viewer_t_scale);
 		folder.Init ( viewer );
		 ph.ToMATLAB ( "temp2.m" );
		 
		ph.FromMATLAB("example3.txt");

		folder.Loop(45);
		if (viewer) SetFloatTime(viewer, viewer_t_scale);
		folder.Init ( viewer );
		
		string tMatlab  = "temp1.m";
		if(!tMatlab.empty()) ph.ToMATLAB ( tMatlab );
					

		delete viewer;
		for ( int i=0; i<nbrCShips; i++ ) delete ship_model[i];
	}

// 	{
// 		lat0 = 54.17057475;
// 		lon0 = 12.10074142;
// 		hgt0 = 0;
// 		
// 		nbrAShips = 0;
// 		nbrAShipsMax = 10;
// 		nbrExtraObstacles = 10;
// 		
// 		withCollReg = 0;
// 		ENDTIME = OPEN_WITH_FIXED_RATIOS;
// 		
// 		viewer_t_scale = 120.0;
// 		
// 		soft_final_constr = 0;
// 		
// 		int nbrCShips = coop_ships_IDs.size();
// 		vector <BaseModel*> ship_model ( nbrCShips );
// 		for(int i=0; i<nbrCShips; i++) ship_model[i] = get_ship_model(coop_ships_IDs[i]);
// 
// 		
// 		GeometricConstr geom_constr(nbrCShips, nbrAShips);
// 
// 		string viewerID = "Optimal Trajectories";
// 		
// 		vector<double> tEnd(tEnd_nu);
// 		vector<double> plow(plow_nu);
// 		vector<double> pupp(pupp_nu);
// 		for(int k=0; k<nbrCShips; k++) {
// 			if(ENDTIME == COMMON_AND_FIXED      ) tEnd[k] = tEnd_nu[k]/scale_t[0];
// 			if(ENDTIME == OPEN_WITH_FIXED_RATIOS) tEnd[k] = tEnd_nu[k]/scale_t[0];
// 			if(ENDTIME == OPEN                  ) tEnd[k] = tEnd_nu[k]/scale_t[k];
// 			plow[k] = plow_nu[k]/scale_t[k];
// 			pupp[k] = pupp_nu[k]/scale_t[k];
// 		}
// 		
// 		vector<int> NmodelS(nbrCShips);
// 		vector<int> NmodelC(nbrCShips);
// 		vector<int> NS(nbrCShips+1);
// 		vector<int> NC(nbrCShips+1);
// 		NS[0] = 0;
// 		NC[0] = 0;
// 
// 		for ( int i=0; i<nbrCShips; i++ ) {
// 			NmodelS[i] = ship_model[i]->nbrS;
// 			NmodelC[i] = ship_model[i]->nbrC;
// 			NS[i+1] = NS[i] + NmodelS[i] + NmodelC[i] + 1;
// 			NC[i+1] = NC[i] + NmodelC[i];
// 		}
// 		
// 		vector <vector<double>> sub_final_nu(nbrCShips);
// 		for(int k=0; k<nbrCShips; k++) sub_final_nu[k].resize(3);
// 
// 		
// 		tw::TWparameter twparameter ( transxml );
// 		map<string,string> args = twparameter.Arguments ( argv,argc );
// 		XMLNode *xml_shipviewer = tw::TWparameter::ReadParams ( shipvxml );
// 		
// 		tw::TWdimension TWdim;
// 		TWdim.ID     = viewerID;
// 		TWdim.n_dis  = twparameter.NDIS;
// 		TWdim.n_ode  = NS[nbrCShips];
// 		TWdim.n_ctrl = NC[nbrCShips];
// 		
// 		if(ENDTIME == COMMON_AND_FIXED      ) TWdim.n_param = 0;
// 		if(ENDTIME == OPEN_WITH_FIXED_RATIOS) TWdim.n_param = 1;
// 		if(ENDTIME == OPEN                  ) TWdim.n_param = nbrCShips;
// 		
// 		TWdim.n_rand  = 0;
// 
// 		
// 		Viewer *viewer = 0;
// 		if ( twparameter.PLOT ) viewer = new tw::Viewer ( &twparameter );
// 		
// 		tw::TWfolder folder ( &twparameter,0 );
// 		
// 		OptTraj ph ( xml_shipviewer, TWdim, nbrCShips, nbrAShips, nbrAShipsMax, nbrObstacles, nbrObstacles_init, nbrObstaclesMax, ship_model, geom_constr, 
// 					distSH, distSO, distSA, distSS, speedbounds,
// 					ENDTIME, tEnd, 
// 					collavoidSH, collavoidSO, collavoidSA, collavoidSS, speedLimit,
// 					Ninterp,  
// 					NS, NC, NmodelS, NmodelC, scale, scale_ctrl, scale_t, viewer_t_scale,
// 					tWgt, eWgt, fWgt, cWgt, plow, pupp, 
// 					ratioToFirst, withCollReg, soft_final_constr, &folder );
// 		ph0 = &ph;
// 		folder0 = &folder;
// 		
// 		ph.setSolver ( &twparameter );
// 		folder.Add ( &ph );
// 		folder.Init();
// 
// 		double x0, y0, z0;
// 		double lat;
// 		double lon;
// 		double psi;
// 		double u;
// 		double v;
// 		double r;
// 		double delta;
// 		double n;
// 		
// 		vector <vector<double>> initial_nu;
// 		vector <vector<double>> final_nu;
// 		initial_nu.push_back({ 54.183615, 12.091236, 160.0, 0, 0, 0, 0, 0, 0 });
// 		final_nu  .push_back({ nan(""), nan(""), 90, nan(""), 0, 0, nan(""), nan(""), nan("") });
// 		
// 		lat   = initial_nu[0][0];
// 		lon   = initial_nu[0][1];
// 		psi   = initial_nu[0][2]*M_PI/180;
// 		u     = initial_nu[0][3];
// 		v     = initial_nu[0][4];
// 		r     = initial_nu[0][5]*M_PI/180;
// 		delta = initial_nu[0][6]*M_PI/180;
// 		n     = initial_nu[0][7]/60;
// 		
// 		LLAtoNED(lat, lon, 0, lat0, lon0, 0, &x0, &y0, &z0);
// 		ph0->var_b_init[0][0] = x0;
// 		ph0->var_b_init[0][1] = y0;
// 		ph0->var_b_init[0][2] = psi;
// 		ph0->var_b_init[0][3] = u;
// 		ph0->var_b_init[0][4] = v;
// 		ph0->var_b_init[0][5] = r;
// 		ph0->var_b_init[0][6] = delta;
// 		ph0->var_b_init[0][7] = n;
// 		for(int i=0; i<8; i++) ph0->var_b_init[0][i] /= ph0->scale[0][i];
// 
// 		ph0->solver->X[ph0->x_index ( 0,ph0->NS[0]+0 )] = x0;
// 		ph0->solver->X[ph0->x_index ( 0,ph0->NS[0]+1 )] = y0;
// 		ph0->solver->X[ph0->x_index ( 0,ph0->NS[0]+2 )] = psi;
// 		ph0->solver->X[ph0->x_index ( 0,ph0->NS[0]+3 )] = u;
// 		ph0->solver->X[ph0->x_index ( 0,ph0->NS[0]+4 )] = v;
// 		ph0->solver->X[ph0->x_index ( 0,ph0->NS[0]+5 )] = r;
// 		ph0->solver->X[ph0->x_index ( 0,ph0->NS[0]+6 )] = delta;
// 		ph0->solver->X[ph0->x_index ( 0,ph0->NS[0]+7 )] = n;
// 		for(int i=0; i<8; i++) ph0->solver->X[ph0->x_index ( 0,ph0->NS[0]+i )] /= ph0->scale[0][i];
// 		
// 		
// 		
// 		lat   = final_nu[0][0];
// 		lon   = final_nu[0][1];
// 		psi   = final_nu[0][2]*M_PI/180;
// 		u     = final_nu[0][3];
// 		v     = final_nu[0][4];
// 		r     = final_nu[0][5]*M_PI/180;
// 		delta = final_nu[0][6]*M_PI/180;
// 		n     = final_nu[0][7]/60;
// 		
// 		LLAtoNED(lat, lon, 0, lat0, lon0, 0, &x0, &y0, &z0);
// 		ph0->var_b_final[0][0] = x0; 
// 		ph0->var_b_final[0][1] = y0; 
// 		ph0->var_b_final[0][2] = psi; 
// 		ph0->var_b_final[0][3] = u; 
// 		ph0->var_b_final[0][4] = v; 
// 		ph0->var_b_final[0][5] = r; 
// 		ph0->var_b_final[0][6] = delta; 
// 		ph0->var_b_final[0][7] = n; 
// 		ph0->var_b_final[0][8] = nan(""); 
// 		for(int i=0; i<9; i++) ph0->var_b_final[0][i] /= ph0->scale[0][i];
// 
// 		
// 		
// 		folder.Init();
// 
// 
// 		ph.FromMATLAB ( "temp1.m" );
// 		
// 	// 	if (viewer) SetFloatTime(viewer, viewer_t_scale);
// 	// 	folder.Init ( viewer );
// 		folder.Loop(0);
// 		if (viewer) SetFloatTime(viewer, viewer_t_scale);
// 		folder.Init ( viewer );
// 		folder.Loop(0);
// 		
// 		string tMatlab  = "temp1.m";
// 		if(!tMatlab.empty()) ph.ToMATLAB ( tMatlab );
// 					
// 
// 		delete viewer;
// 		delete xml_shipviewer;
// 		for ( int i=0; i<nbrCShips; i++ ) delete ship_model[i];
// 	}
	
	return 0;
}



