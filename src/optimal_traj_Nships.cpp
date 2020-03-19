#include "optimal_traj_Nships.h"
#include "global.h"


OptTraj::OptTraj ( 
	XMLNode *xmlmain,
	const tw::TWdimension &TWdata,
	int nbrCShips,
	int nbrAShips,
	int nbrAShipsMax,
	int nbrObstacles,
	int nbrObstacles_init,
	int nbrObstaclesMax,
	const vector <BaseModel*> &ship_model,
	const GeometricConstr &geom_constr,
	vector<double>         distSH,
	vector<double>         distSO,
	vector<vector<double>> distSA,
	vector<vector<double>> distSS,
	vector<double> speedbounds,
	EndTimeBool ENDTIME,
	vector<double> tEnd,
	bool collavoidSH,
	bool collavoidSO,
	bool collavoidSA,
	bool collavoidSS,
	bool speedLimit,
	int Ninterp,
	vector<int> NS,
	vector<int> NC,
	vector<int> NmodelS,
	vector<int> NmodelC,
	vector <vector<double>> scale,
	vector <vector<double>> scale_ctrl,
	vector<double> scale_t, 
	double viewer_t_scale,
	vector<double> tWgt,
	vector<double> eWgt,
	vector<vector<double>> fWgt,
	vector<vector<double>> cWgt,
	vector<double> plow,
	vector<double> pupp,
	const vector<double> ratioToFirst,
	bool withCollReg,
	bool soft_final_constr,
	tw::TWfolder* folder,
	string resultfile

	)
	:
	TransWorhpProblem ( TWdata ) ,
	nbrCShips(nbrCShips) , 
	nbrAShips(nbrAShips) , 
	nbrAShipsMax(nbrAShipsMax),
	nbrObstacles(nbrObstacles) , 
	nbrObstacles_init(nbrObstacles_init),
	nbrObstaclesMax(nbrObstaclesMax),
	ship_model ( ship_model ) ,
	geom_constr ( geom_constr ),
	distSH(distSH),
	distSO(distSO),
	distSA(distSA),
	distSS(distSS),
	speedbounds(speedbounds),
	ENDTIME (ENDTIME), 
	tEnd (tEnd), 
	collavoidSH ( collavoidSH ),
	collavoidSO ( collavoidSO ),
	collavoidSA ( collavoidSA ),
	collavoidSS ( collavoidSS ),
	speedLimit (speedLimit),
	Ninterp (Ninterp),
	NS(NS),
	NC(NC),
	NmodelS(NmodelS),
	NmodelC(NmodelC),
	scale(scale),
	scale_ctrl(scale_ctrl),
	scale_t(scale_t),
	viewer_t_scale(viewer_t_scale),
	tWgt(tWgt),
	eWgt(eWgt),
	fWgt(fWgt),
	cWgt(cWgt),
	plow(plow),
	pupp(pupp),
	ratioToFirst( ratioToFirst ),
	withCollReg (withCollReg),
	soft_final_constr(soft_final_constr),
	folder_ptr(folder),
	resultfile(resultfile) {

	if ( xmlmain ) scenenode = xmlmain->GetFirstChild ( "SCENE" );
	
	ctrl   =  new double[ NC[nbrCShips] ];
	state  =  new double[ NS[nbrCShips] ];
	dstate =  new double[ NS[nbrCShips] ];
	
	var_b_init.resize(nbrCShips);
	var_b_final.resize(nbrCShips);
	final_s.resize(nbrCShips);
	for(int k=0; k<nbrCShips; k++) {
		var_b_init[k].resize(NS[k+1]-NS[k]);
		for(int i=0; i<NS[k+1]-NS[k]; i++) var_b_init[k][i] = 0;
	}
	for(int k=0; k<nbrCShips; k++) {
		var_b_final[k].resize(NS[k+1]-NS[k]);
		for(int i=0; i<NS[k+1]-NS[k]; i++) var_b_final[k][i] = 0;
	}
	for(int k=0; k<nbrCShips; k++) {
		final_s[k].resize(NS[k+1]-NS[k]);
		for(int i=0; i<NS[k+1]-NS[k]; i++) final_s[k][i] = 0;
	}
	
	AISShips.resize(nbrAShips);	
}

void OptTraj::OpenWindows ( tw::Viewer *viewer ) {
	viewer->ThreeD ( "Ship tw::Viewer", scenenode, shipviewer3d );
}

void OptTraj::selectWindows ( tw::Viewer *viewer ) {
	int idx = 0;
	for ( int k=0; k<nbrCShips; k++ ) {
		for ( int i=0; i<NmodelS[k]; i++ ) {
// 					viewer->AddStateView ( NS[k]+i , "" );
			viewer->Data("", scaleplot, idx);
			idx++;
		}
		for ( int i=0; i<NmodelC[k]; i++ ) {
// 					viewer->AddStateView ( NS[k]+NmodelS[k]+i , "" );
			viewer->Data("", scaleplot, idx);
			idx++;
		}
// 				viewer->AddStateView ( NS[k]+NmodelS[k]+NmodelC[k] , "" );
		viewer->Data("", scaleplot, idx);
		idx++;
	}
// 			viewer->AddStateView ( NS[nbrCShips], "" );
// 			viewer->Data("", scaleplot, idx);
// 			idx++;
	for ( int k=0; k<nbrCShips; k++ ) {
		for ( int i=0; i<NmodelC[k]; i++ ) {
// 					viewer->AddControlView ( NC[k]+i , "" );
			viewer->Data("", scaleplot, idx);
			idx++;
		}
	}
}

void OptTraj::p_init ( double *p ) {
	if(ENDTIME == COMMON_AND_FIXED      ) /*nothing*/                                  ;
	if(ENDTIME == OPEN_WITH_FIXED_RATIOS) p[0] = tEnd[0]                               ;
	if(ENDTIME == OPEN                  ) for(int k=0; k<nbrCShips; k++) p[k] = tEnd[k];
}

double OptTraj::obj() {
	double objective = 0;
// 	double obj_prev = 0;
	if(ENDTIME == COMMON_AND_FIXED      ) /*nothing*/      ;
	if(ENDTIME == OPEN_WITH_FIXED_RATIOS) objective += tWgt[0]*p(0);
	
// 	cout << "obj: " << objective - obj_prev << endl;
// 	obj_prev = objective;
	
	for ( int k=0; k<nbrCShips; k++ ) {
		if(ENDTIME == OPEN) objective += tWgt[k]*p(k);
		objective += eWgt[k]*x( n_dis-1, NS[k]+NmodelS[k]+NmodelC[k] );
	}
	
// 	cout << "obj: " << objective - obj_prev << endl;
// 	obj_prev = objective;
	
	if(soft_final_constr) {
		for ( int k=0; k<nbrCShips; k++ ) {
			for ( int i=0; i<NmodelS[k]+NmodelC[k]; i++ ) {
				if ( !std::isnan(fWgt[k][i]) && !std::isnan(var_b_final[k][i]) ) {
					objective += fWgt[k][i]*(x(n_dis-1,NS[k]+i) - var_b_final[k][i]) * (x(n_dis-1,NS[k]+i) - var_b_final[k][i]);
				}
			}
		}
	}
	
// 	cout << "obj: " << objective - obj_prev << endl;
// 	cout << endl;
// 	obj_prev = objective;

	objective /= n_dis;

	return objective;
}

bool OptTraj::obj_structure ( tw::DiffStructure &s ) {

	if(ENDTIME == COMMON_AND_FIXED      ) /*nothing*/          ;
	if(ENDTIME == OPEN_WITH_FIXED_RATIOS) s( 0, p_index ( 0 ) );
	
	for ( int k=0; k<nbrCShips; k++ ) {
		if(ENDTIME == OPEN) s( 0, p_index ( k ) );
		s ( 0, x_index( n_dis-1, NS[k]+NmodelS[k]+NmodelC[k] ) );
	}
	
	for ( int k=0; k<nbrCShips; k++ ) {
		for ( int i=0; i<NmodelS[k]+NmodelC[k]; i++ ) {
			if ( !std::isnan(fWgt[k][i]) && !std::isnan(var_b_final[k][i]) ) {
				s ( 0, x_index( n_dis-1, NS[k]+i ) );
			}
		}
	}
	
	return true;
}

bool OptTraj::obj_diff ( tw::DiffStructure &s ) {
	return false;
}

void OptTraj::ode ( double *dx, double t, const double *x, const double *u, const double *p ) {
	for ( int k=0; k<nbrCShips; k++ ) {
		double multiplier;
		if(ENDTIME == COMMON_AND_FIXED      ) multiplier =                 tEnd[0]*scale_t[0];
		if(ENDTIME == OPEN_WITH_FIXED_RATIOS) multiplier = ratioToFirst[k]*p[0]   *scale_t[0];
		if(ENDTIME == OPEN)                   multiplier =                 p[k]   *scale_t[k];

		for ( int i=0; i<NmodelS[k]; i++ ) state[i] = x[NS[k]           +i]*scale[k][i];
		for ( int i=0; i<NmodelC[k]; i++ ) ctrl [i] = x[NS[k]+NmodelS[k]+i]*scale[k][NmodelS[k]+i];

		double q[9] = {0,0,0,0,0,0,0,0,1000.};

		ship_model[k]->ode(dstate, state, ctrl, NULL, q);

		for ( int i=0; i<NmodelS[k]; i++ ) dx[NS[k]           +i] = dstate[i]  * multiplier / scale[k][i];
		for ( int i=0; i<NmodelC[k]; i++ ) dx[NS[k]+NmodelS[k]+i] = scale_ctrl[k][i]*u[NC[k]+i] * multiplier / scale[k][NmodelS[k]+i];

		dx[NS[k]+NmodelS[k]+NmodelC[k]]  = 0;
		dx[NS[k]+NmodelS[k]+NmodelC[k]] += cWgt[k][0]*scale_ctrl[k][0]*u[NC[k]+0]*scale_ctrl[k][0]*u[NC[k]+0] * multiplier / scale[k][NmodelS[k]+NmodelC[k]];
		dx[NS[k]+NmodelS[k]+NmodelC[k]] += cWgt[k][1]*scale_ctrl[k][1]*u[NC[k]+1]*scale_ctrl[k][1]*u[NC[k]+1] * multiplier / scale[k][NmodelS[k]+NmodelC[k]];
		dx[NS[k]+NmodelS[k]+NmodelC[k]] += cWgt[k][2]*x[NS[k]+7]*x[NS[k]+7] * multiplier / scale[k][NmodelS[k]+NmodelC[k]];
	}
}

bool OptTraj::ode_structure ( tw::DiffStructure &s ) {
	
	for ( int k=0; k<nbrCShips; k++ ) {
		ship_model[k]->ode_structure ( s, NS[k] );
		if(ENDTIME == COMMON_AND_FIXED      ) /*nothing*/                                                      ;
		if(ENDTIME == OPEN_WITH_FIXED_RATIOS) for ( int i=0; i<NmodelS[k]; i++ ) s( NS[k]+i, p_indexode ( 0 ) );
		if(ENDTIME == OPEN)                   for ( int i=0; i<NmodelS[k]; i++ ) s( NS[k]+i, p_indexode ( k ) );

		for ( int i=0; i<NmodelC[k]; i++ ) {
			s ( NS[k]+NmodelS[k]+i, u_indexode ( NC[k]+i ) );
			if(ENDTIME == COMMON_AND_FIXED      ) /*nothing*/                               ;
			if(ENDTIME == OPEN_WITH_FIXED_RATIOS) s ( NS[k]+NmodelS[k]+i, p_indexode ( 0 ) );
			if(ENDTIME == OPEN)                   s ( NS[k]+NmodelS[k]+i, p_indexode ( k ) );
			s ( NS[k]+NmodelS[k]+NmodelC[k], u_indexode ( NC[k]+i ) );
			s ( NS[k]+NmodelS[k]+NmodelC[k], x_indexode ( NS[k]+7 ) );
		}
		if(ENDTIME == COMMON_AND_FIXED      ) /*nothing*/                                        ;
		if(ENDTIME == OPEN_WITH_FIXED_RATIOS) s ( NS[k]+NmodelS[k]+NmodelC[k], p_indexode ( 0 ) );
		if(ENDTIME == OPEN)                   s ( NS[k]+NmodelS[k]+NmodelC[k], p_indexode ( k ) );
	}
	return true;
}

void OptTraj::p_boundary ( double *p_low, double *p_upp ) {
	if(ENDTIME == COMMON_AND_FIXED      ) /*nothing*/;
	if(ENDTIME == OPEN_WITH_FIXED_RATIOS) {
		p_low[0] = plow[0];
		p_upp[0] = pupp[0];
	}
	if(ENDTIME == OPEN) {
		for(int k=0; k<nbrCShips; k++) {
			p_low[k] = plow[k];
			p_upp[k] = pupp[k];
		}
	}
}

void OptTraj::x_boundary ( double *x_low, double *x_upp ) {
	for ( int k=0; k<nbrCShips; k++ ) {
		double *low = new double[NmodelS[k]+NmodelC[k]];
		double *upp = new double[NmodelS[k]+NmodelC[k]];
		ship_model[k]->get_bounds ( low, upp );

		for ( int i=0; i<NmodelS[k]+NmodelC[k]; i++ ) {
			if ( !std::isnan ( low[i] ) ) x_low[NS[k]+i] = low[i]/scale[k][i];

			if ( !std::isnan ( upp[i] ) ) x_upp[NS[k]+i] = upp[i]/scale[k][i];
		}
		delete[] low;
		delete[] upp;
	}
}

void OptTraj::u_boundary ( double *u_low, double *u_upp ) {
	for ( int k=0; k<nbrCShips; k++ ) {
		double *low = new double[NmodelC[k]];
		double *upp = new double[NmodelC[k]];
		ship_model[k]->get_diff_bounds ( low, upp );

		for ( int i=0; i<NmodelC[k]; i++ ) {
			if ( !std::isnan ( low[i] ) ) u_low[NC[k]+i] = low[i]/scale_ctrl[k][i];

			if ( !std::isnan ( upp[i] ) ) u_upp[NC[k]+i] = upp[i]/scale_ctrl[k][i];
		}
		delete[] low;
		delete[] upp;
	}
}

void OptTraj::var_boundary ( double *x_low, double *x_upp ) {
	for ( int k=0; k<nbrCShips; k++ ) {
		// interface conditions for the creation of motion primitives 
				
				/********************** 
				x_low[u_index(0,0)] = x_upp[u_index(0,0)] = 0.0;
				x_low[u_index(n_dis-1,0)] = x_upp[u_index(n_dis-1,0)] = 0.0;
				
				x_low[u_index(0,1)] = x_upp[u_index(0,1)] = 0.0;
				x_low[u_index(n_dis-1,1)] = x_upp[u_index(n_dis-1,1)] = 0.0;
				*************************/
		for ( int i=0; i<NmodelS[k]+NmodelC[k]+1; i++ ) {
			if ( !std::isnan ( var_b_init[k][i] ) ) {
				x_low[x_index ( 0,NS[k]+i ) ] = x_upp[x_index ( 0,NS[k]+i ) ] = var_b_init[k][i];
			}
		}
	}
	
	if(!soft_final_constr) {
		for ( int k=0; k<nbrCShips; k++ ) {
			for ( int i=0; i<NmodelS[k]+NmodelC[k]; i++ ) {
				if ( !std::isnan ( var_b_final[k][i] ) ) {
					//cout << "§§§  " << i << "  " << scale[k][i]*var_b_final[k][i] << endl;
					x_low[x_index ( n_dis-1,NS[k]+i ) ] = x_upp[x_index ( n_dis-1,NS[k]+i ) ] = var_b_final[k][i];
				}
			}
		}
	}
	

		for ( int k=0; k<nbrCShips; k++ ) {
			for ( int i=0; i<NmodelS[k]+NmodelC[k]; i++ ) {

					//x_low[x_index ( n_dis-1,NS[k]+3 ) ] = x_upp[x_index ( n_dis-1,NS[k]+3 ) ] = 0.01;

			}
		}

}
		
		
void OptTraj::rand(double *r) {
	double sqdist, sqarea;
	double ship_area;
	double x0, y0, psi;
	double dt, dtfull;
	int idx=0;
	
	if(collavoidSH) {
		for(int k=0; k<nbrCShips; k++) {
			for(int i=0; i<n_dis-1; i++) {
				for(int j=0; j<Ninterp+1; j++) {
					dt = (double)j/(Ninterp+1);
					if(ENDTIME == COMMON_AND_FIXED      ) dtfull = tEnd[0]/(n_dis-1)*scale_t[0];
					if(ENDTIME == OPEN_WITH_FIXED_RATIOS) dtfull = p(0)*ratioToFirst[k]/(n_dis-1)*scale_t[0];
					if(ENDTIME == OPEN)                   dtfull = p(k)/(n_dis-1)*scale_t[k];
					cubicSpline(i, k, dt, dtfull, 0, &x0, &y0, &psi);
					x0  *= scale[k][0];
					y0  *= scale[k][1];
					psi *= scale[k][2];
					geom_constr.updateCShip(k, x0, y0, psi);
					geom_constr.ConstrCShipHarbr( k, 0, &sqdist, &sqarea, geom_constr.harb_polygons[0].map_grid.active );
					ship_area = geom_constr.getAreaOfShip(k);
					r[idx] = 1e-4*(-sqdist + sqarea - geom_constr.getAFactor()*ship_area*ship_area);
					idx++;
				}
			}
			dt = 0; 
			dtfull = 0;
			cubicSpline(n_dis-1, k, dt, dtfull, 0, &x0, &y0, &psi);
			x0  *= scale[k][0];
			y0  *= scale[k][1];
			psi *= scale[k][2];
			geom_constr.updateCShip(k, x0, y0, psi);
			geom_constr.ConstrCShipHarbr( k, 0, &sqdist, &sqarea, geom_constr.harb_polygons[0].map_grid.active );
			ship_area = geom_constr.getAreaOfShip(k);
			r[idx] = 1e-4*(-sqdist + sqarea - geom_constr.getAFactor()*ship_area*ship_area);
			idx++;
		}
	}
	if(collavoidSO) {
		for(int k=0; k<nbrCShips; k++) {
			for ( int l=0; l<nbrObstaclesMax; l++ ) {
				if(l<nbrObstacles) {
					for(int i=0; i<n_dis-1; i++) {
						for(int j=0; j<Ninterp+1; j++) {
							dt = (double)j/(Ninterp+1);
							if(ENDTIME == COMMON_AND_FIXED      ) dtfull = tEnd[0]/(n_dis-1)*scale_t[0];
							if(ENDTIME == OPEN_WITH_FIXED_RATIOS) dtfull = p(0)*ratioToFirst[k]/(n_dis-1)*scale_t[0];
							if(ENDTIME == OPEN)                   dtfull = p(k)/(n_dis-1)*scale_t[k];
							cubicSpline(i, k, dt, dtfull, 0, &x0, &y0, &psi);
							x0  *= scale[k][0];
							y0  *= scale[k][1];
							psi *= scale[k][2];
							geom_constr.updateCShip(k, x0, y0, psi);
							geom_constr.ConstrCShipObsta( k, l, &sqdist, &sqarea, geom_constr.obst_polygons[l].map_grid.active);
							r[idx] = 1e-0*(sqdist - sqarea);
							idx++;
						}
					}
					dt = 0;
					dtfull = 0;
					cubicSpline(n_dis-1, k, dt, dtfull, 0, &x0, &y0, &psi);
					x0  *= scale[k][0];
					y0  *= scale[k][1];
					psi *= scale[k][2];
					geom_constr.updateCShip(k, x0, y0, psi);
					geom_constr.ConstrCShipObsta( k, l, &sqdist, &sqarea, geom_constr.obst_polygons[l].map_grid.active);
					r[idx] = 1e-0*(sqdist - sqarea);
					idx++;
				} else {
					r[idx] = 0;
					idx++;
				}
			}
		}
	}
	if(collavoidSA) {
		for(int k=0; k<nbrCShips; k++) {
			for ( int l=0; l<nbrAShipsMax; l++ ) {
				if(l<nbrAShips) {
					for(int i=0; i<n_dis-1; i++) {
						for(int j=0; j<Ninterp+1; j++) {
							dt = (double)j/(Ninterp+1);
							if(ENDTIME == COMMON_AND_FIXED      ) dtfull = tEnd[0]/(n_dis-1)*scale_t[0];
							if(ENDTIME == OPEN_WITH_FIXED_RATIOS) dtfull = p(0)*ratioToFirst[k]/(n_dis-1)*scale_t[0];
							if(ENDTIME == OPEN)                   dtfull = p(k)/(n_dis-1)*scale_t[k];
							cubicSpline(i, k, dt, dtfull, 0, &x0, &y0, &psi);
							x0  *= scale[k][0];
							y0  *= scale[k][1];
							psi *= scale[k][2];
							geom_constr.updateCShip(k, x0, y0, psi);
							if(ENDTIME == COMMON_AND_FIXED      ) dt = tEnd[0]*scale_t[0]*i/(n_dis-1)+tEnd[0]*scale_t[0]/(n_dis-1)*(j+1.)/(Ninterp+1);
							if(ENDTIME == OPEN_WITH_FIXED_RATIOS) dt = ratioToFirst[k]*p(0)*scale_t[0]*i/(n_dis-1) + ratioToFirst[k]*p(0)*scale_t[0]/(n_dis-1)*(j+1.)/(Ninterp+1);
							if(ENDTIME == OPEN)                   dt = p(k)*scale_t[k]*i/(n_dis-1)+p(k)*scale_t[k]/(n_dis-1)*(j+1.)/(Ninterp+1);

							x0  = AISShips[l].x0 + dt*AISShips[l].sog*cos(AISShips[l].cog);
							y0  = AISShips[l].y0 + dt*AISShips[l].sog*sin(AISShips[l].cog);
							psi = AISShips[l].psi;
							geom_constr.updateAShip(l, x0, y0, psi);

							geom_constr.ConstrCShipAShip( k, l, &sqdist, &sqarea, 0);
							r[idx] = 1e-0*(sqdist - sqarea);
							idx++;
						}
					}
					dt = 0;
					dtfull = 0;
					cubicSpline(n_dis-1, k, dt, dtfull, 0, &x0, &y0, &psi);
					x0  *= scale[k][0];
					y0  *= scale[k][1];
					psi *= scale[k][2];
					geom_constr.updateCShip(k, x0, y0, psi);
					if(ENDTIME == COMMON_AND_FIXED      ) dt = tEnd[0]*scale_t[0];
					if(ENDTIME == OPEN_WITH_FIXED_RATIOS) dt = ratioToFirst[k]*p(0)*scale_t[0];
					if(ENDTIME == OPEN)                   dt = p(k)*scale_t[k];

					x0  = AISShips[l].x0 + dt*AISShips[l].sog*cos(AISShips[l].cog);
					y0  = AISShips[l].y0 + dt*AISShips[l].sog*sin(AISShips[l].cog);
					psi = AISShips[l].psi;
					geom_constr.updateAShip(l, x0, y0, psi);

					geom_constr.ConstrCShipAShip( k, l, &sqdist, &sqarea, 0);
					r[idx] = 1e-0*(sqdist - sqarea);
					idx++;
				} else {
					r[idx] = 0;
					idx++;
				}
			}
		}
	}
	
	if(collavoidSS && (ENDTIME == COMMON_AND_FIXED || ENDTIME == OPEN_WITH_FIXED_RATIOS)) {
		double t1, dt1, dt2;
		int ilow;
		int l,k;
		for ( int kk=0; kk<nbrCShips; kk++ ) {
			for ( int ll=kk+1; ll<nbrCShips; ll++ ) {
				if(ENDTIME == COMMON_AND_FIXED      ) {l = ll; k = kk;}
				if(ENDTIME == OPEN_WITH_FIXED_RATIOS) {
					if(ratioToFirst[kk] > ratioToFirst[ll]) {
						l = kk; k = ll;
					} else { 
						l = ll; k = kk;
					}
				}
				if(ENDTIME == OPEN) {
					if(p(kk)*scale_t[kk] > p(ll)*scale_t[ll]) {
						l = kk; k = ll;
					} else { 
						l = ll; k = kk;
					}
				}
				for(int i=0; i<n_dis-1; i++) {
					for(int j=0; j<Ninterp+1; j++) {
						dt = (double)j/(Ninterp+1);
						if(ENDTIME == COMMON_AND_FIXED      ) dtfull = tEnd[0]*scale_t[0]/(n_dis-1);
						if(ENDTIME == OPEN_WITH_FIXED_RATIOS) dtfull = p(0)*scale_t[0]*ratioToFirst[l]/(n_dis-1);
						cubicSpline(i, l, dt, dtfull, 0, &x0, &y0, &psi);
						x0  *= scale[k][0];
						y0  *= scale[k][1];
						psi *= scale[k][2];
						geom_constr.updateCShip(l, x0, y0, psi);
						
						if(ENDTIME == COMMON_AND_FIXED      ) {dt1 = tEnd[0]*scale_t[0]/(n_dis-1); dt2 = tEnd[0]*scale_t[0]/(n_dis-1);}
						if(ENDTIME == OPEN_WITH_FIXED_RATIOS) {dt1 = ratioToFirst[l]*p(0)*scale_t[0]/(n_dis-1); dt2 = ratioToFirst[k]*p(0)*scale_t[0]/(n_dis-1);}
						
						t1 = dt1*(i+1.*j/(Ninterp+1));
						ilow = (int) floor(t1/dt2);
						
						dt = t1 - dt2*ilow;
						
						if(ilow+1 <= n_dis-1) {
							cubicSpline(ilow, k, dt/dt2, dt2, 0, &x0, &y0, &psi);
							x0  *= scale[k][0];
							y0  *= scale[k][1];
							psi *= scale[k][2];
							geom_constr.updateCShip(k, x0, y0, psi);
							geom_constr.ConstrCShipCShip( k, l, &sqdist, &sqarea, 0);
							r[idx] = 1e-0*(sqdist - sqarea);
							idx++;
						} else {
							if(ENDTIME == COMMON_AND_FIXED      ) { dt = t1 - tEnd[0]*scale_t[0];              dtfull = tEnd[0]*scale_t[0];             }
							if(ENDTIME == OPEN_WITH_FIXED_RATIOS) { dt = t1 - ratioToFirst[k]*p(0)*scale_t[0]; dtfull = ratioToFirst[k]*p(0)*scale_t[0];}
							cubicSpline(n_dis-1, k, dt, dtfull, 1, &x0, &y0, &psi);
							x0  *= scale[k][0];
							y0  *= scale[k][1];
							psi *= scale[k][2];
							geom_constr.updateCShip(k, x0, y0, psi);
							geom_constr.ConstrCShipCShip( k, l, &sqdist, &sqarea, 0);
							r[idx] = 1e-0*(sqdist - sqarea);
							idx++;
						}
					}
				}
				dt = 0;
				dtfull = 0;
				cubicSpline(n_dis-1, l, dt, dtfull, 0, &x0, &y0, &psi);
				x0  *= scale[k][0];
				y0  *= scale[k][1];
				psi *= scale[k][2];
				geom_constr.updateCShip(l, x0, y0, psi);
				
				if(ENDTIME == COMMON_AND_FIXED      ) {dt1 = tEnd[0]*scale_t[0]/(n_dis-1); dt2 = tEnd[0]*scale_t[0]/(n_dis-1);                          }
				if(ENDTIME == OPEN_WITH_FIXED_RATIOS) {dt1 = ratioToFirst[l]*p(0)*scale_t[0]/(n_dis-1); dt2 = ratioToFirst[k]*p(0)*scale_t[0]/(n_dis-1);}
				
				t1 = dt1*(n_dis-1);
				ilow = (int) floor(t1/dt2);
				
				dt = t1 - dt2*ilow;
				
				if(ilow+1 <= n_dis-1) {
					cubicSpline(ilow, k, dt/dt2, dt2, 0, &x0, &y0, &psi);
					x0  *= scale[k][0];
					y0  *= scale[k][1];
					psi *= scale[k][2];
					geom_constr.updateCShip(k, x0, y0, psi);
					geom_constr.ConstrCShipCShip( k, l, &sqdist, &sqarea, 0);
					r[idx] = 1e-0*(sqdist - sqarea);
					idx++;
				} else {
					if(ENDTIME == COMMON_AND_FIXED      ) { dt = t1 - tEnd[0]*scale_t[0];              dtfull = tEnd[0]*scale_t[0];             }
					if(ENDTIME == OPEN_WITH_FIXED_RATIOS) { dt = t1 - ratioToFirst[k]*p(0)*scale_t[0]; dtfull = ratioToFirst[k]*p(0)*scale_t[0];}
					cubicSpline(n_dis-1, k, dt, dtfull, 1, &x0, &y0, &psi);
					x0  *= scale[k][0];
					y0  *= scale[k][1];
					psi *= scale[k][2];
					geom_constr.updateCShip(k, x0, y0, psi);
					geom_constr.ConstrCShipCShip( k, l, &sqdist, &sqarea, 0);
					r[idx] = 1e-0*(sqdist - sqarea);
					idx++;
				}
			}
		}
	}
	if(speedLimit) {
		for(int k=0; k<nbrCShips; k++) {
				for(int i=0; i<n_dis; i++) {
					x0  = x(i, NS[k]+0)*scale[k][0];
					y0  = x(i, NS[k]+1)*scale[k][1];
					psi = x(i, NS[k]+2)*scale[k][2];
					geom_constr.updateCShip(k, x0, y0, psi);
					r[idx] = 1e-0*( geom_constr.SpeedBound(k, speedbounds) - x(i, NS[k]+3)*scale[k][3] );
					idx++;
				}
		}
	}
}
		
void OptTraj::rand_boundary(double *r_low, double *r_upp) {
	
	int idx=0;
	
	if(collavoidSH) {
		for(int k=0; k<nbrCShips; k++) {
			for(int i=0; i<n_dis-1; i++) {
				for(int j=0; j<Ninterp+1; j++) {
					r_low[idx]=1e-4*distSH[k]*distSH[k];
// 					r_low[idx]=1e-2*distSH[k]*tanh(distSH[k]);
					r_upp[idx]=1e20;
					idx++;
				}
			}
			r_low[idx]=1e-4*distSH[k]*distSH[k];
// 			r_low[idx]=1e-2*distSH[k]*tanh(distSH[k]);
			r_upp[idx]=1e20;
			idx++;					
		}
	}
	if(collavoidSO) {
		for(int k=0; k<nbrCShips; k++) {
			for ( int l=0; l<nbrObstaclesMax; l++ ) {
				if(l<nbrObstacles) {
					for(int i=0; i<n_dis-1; i++) {
						for(int j=0; j<Ninterp+1; j++) {
							r_low[idx]=1e-0*distSO[k]*distSO[k];
							r_upp[idx]=1e20;
							idx++;
						}
					}
					r_low[idx]=1e-0*distSO[k]*distSO[k];
					r_upp[idx]=1e20;
					idx++;
				} else {
					r_low[idx]=0;
					r_upp[idx]=0;
					idx++;
				}
			}
		}
	}
	if(collavoidSA) {
		for(int k=0; k<nbrCShips; k++) {
			for ( int l=0; l<nbrAShipsMax; l++ ) {
				if(l<nbrAShips) {
					for(int i=0; i<n_dis-1; i++) {
						for(int j=0; j<Ninterp+1; j++) {
							r_low[idx]=1e-0*distSA[k][l]*distSA[k][l];
							r_upp[idx]=1e20;
							idx++;
						}
					}
					r_low[idx]=1e-0*distSA[k][l]*distSA[k][l];
					r_upp[idx]=1e20;
					idx++;
				} else {
					r_low[idx]=0;
					r_upp[idx]=0;
					idx++;
				}
			}
		}
	}
	if(collavoidSS && (ENDTIME == COMMON_AND_FIXED || ENDTIME == OPEN_WITH_FIXED_RATIOS)) {
		for ( int k=0; k<nbrCShips; k++ ) {
			for ( int l=k+1; l<nbrCShips; l++ ) {
				for(int i=0; i<n_dis-1; i++) {
					for(int j=0; j<Ninterp+1; j++) {
						r_low[idx]=1e-0*distSS[k][l-(k+1)]*distSS[k][l-(k+1)];
						r_upp[idx]=1e20;
						idx++;
					}
				}
				r_low[idx]=1e-0*distSS[k][l-(k+1)]*distSS[k][l-(k+1)];
				r_upp[idx]=1e20;
				idx++;
			}
		}
	}
	if(speedLimit) {
		for(int k=0; k<nbrCShips; k++) {
			for ( int l=0; l<1; l++ ) {
				for(int i=0; i<n_dis; i++) {
					r_low[idx]=0;
					r_upp[idx]=1e20;
					idx++;
				}
			}
		}
	}
}
		
bool OptTraj::rand_structure(tw::DiffStructure &s) {
	int idx=0;
	if(collavoidSH) {
		for(int k=0; k<nbrCShips; k++) {
			for(int i=0; i<n_dis-1; i++) {
				for(int j=0; j<Ninterp+1; j++) {
					if(j!=0) {
						cubicSpline_structure(s, i, k, idx, NORMAL);
						if(ENDTIME == COMMON_AND_FIXED      ) /*nothing*/;
						if(ENDTIME == OPEN_WITH_FIXED_RATIOS) s ( idx, p_index(0) );
						if(ENDTIME == OPEN)                   s ( idx, p_index(k) );
					}
					else {
						cubicSpline_structure(s, i, k, idx, SINGLEPOINT);
					}
					idx++;
				}
			}
			cubicSpline_structure(s, n_dis-1, k, idx, SINGLEPOINT);
			idx++;
		}
	}
	if(collavoidSO) {
		for(int k=0; k<nbrCShips; k++) {
			for ( int l=0; l<nbrObstaclesMax; l++ ) {
				if(l<nbrObstacles) {
					for(int i=0; i<n_dis-1; i++) {
						for(int j=0; j<Ninterp+1; j++) {
							if(j!=0) {
								cubicSpline_structure(s, i, k, idx, NORMAL);
								if(ENDTIME == COMMON_AND_FIXED      ) /*nothing*/;
								if(ENDTIME == OPEN_WITH_FIXED_RATIOS) s ( idx, p_index(0) );
								if(ENDTIME == OPEN)                   s ( idx, p_index(k) );
							}
							else {
								cubicSpline_structure(s, i, k, idx, SINGLEPOINT);
							}
							idx++;
						}
					}
					cubicSpline_structure(s, n_dis-1, k, idx, SINGLEPOINT);
					idx++;
				} else {
					idx++;
				}
			}
		}
	}
	if(collavoidSA) {
		for(int k=0; k<nbrCShips; k++) {
			for ( int l=0; l<nbrAShipsMax; l++ ) {
				if(l<nbrAShips) {
					for(int i=0; i<n_dis-1; i++) {
						for(int j=0; j<Ninterp+1; j++) {
							if(j!=0) {
								cubicSpline_structure(s, i, k, idx, NORMAL);
							}
							else {
								cubicSpline_structure(s, i, k, idx, SINGLEPOINT);
							}
							if(ENDTIME == COMMON_AND_FIXED      ) /*nothing*/             ;
							if(ENDTIME == OPEN_WITH_FIXED_RATIOS) s ( idx, p_index ( 0 ) );
							if(ENDTIME == OPEN)                   s ( idx, p_index ( k ) );
							idx++;
						}
					}
					cubicSpline_structure(s, n_dis-1, k, idx, SINGLEPOINT);
					if(ENDTIME == COMMON_AND_FIXED      ) /*nothing*/             ;
					if(ENDTIME == OPEN_WITH_FIXED_RATIOS) s ( idx, p_index ( 0 ) );
					if(ENDTIME == OPEN)                   s ( idx, p_index ( k ) );
					idx++;
				} else {
					idx++;
				}
			}
		}
	}
	if(collavoidSS && (ENDTIME == COMMON_AND_FIXED || ENDTIME == OPEN_WITH_FIXED_RATIOS)) {
		double t1, dt1, dt2, dt;
		int ilow;
		int k,l;
		for ( int kk=0; kk<nbrCShips; kk++ ) {
			for ( int ll=kk+1; ll<nbrCShips; ll++ ) {
				if(ENDTIME == COMMON_AND_FIXED      ) {l = ll; k = kk;}
				if(ENDTIME == OPEN_WITH_FIXED_RATIOS) {
					if(ratioToFirst[kk] > ratioToFirst[ll]) {
						l = kk; k = ll;
					} else { 
						l = ll; k = kk;
					}
				}

				for(int i=0; i<n_dis-1; i++) {
					for(int j=0; j<Ninterp+1; j++) {

						if(ENDTIME == COMMON_AND_FIXED      ) ilow = (int) floor((i+1.*j/(Ninterp+1)));
						if(ENDTIME == OPEN_WITH_FIXED_RATIOS) ilow = (int) floor((i+1.*j/(Ninterp+1))*ratioToFirst[l]/ratioToFirst[k]);
						
						if(j!=0) {
							cubicSpline_structure(s, i, l, idx, NORMAL);
							if(ENDTIME == OPEN_WITH_FIXED_RATIOS) s ( idx, p_index ( 0 ) );
						}
						else {
							cubicSpline_structure(s, i, l, idx, SINGLEPOINT);
						}

						if(ilow+1 <= n_dis-1) {
							if(i!=0 || j!=0) {
								cubicSpline_structure(s, ilow, k, idx, NORMAL);
								if(ENDTIME == OPEN_WITH_FIXED_RATIOS) s ( idx, p_index ( 0 ) );
							} else {
								cubicSpline_structure(s, ilow, k, idx, SINGLEPOINT);
							}
							idx++;
						} else {
							cubicSpline_structure(s, n_dis-1, k, idx, EXTRAPOLATE);
							if(ENDTIME == OPEN_WITH_FIXED_RATIOS) s ( idx, p_index ( 0 ) );
							idx++;
						}
					}
				}
				if(ENDTIME == COMMON_AND_FIXED      ) ilow = (int) floor(n_dis);
				if(ENDTIME == OPEN_WITH_FIXED_RATIOS) ilow = (int) floor(n_dis*ratioToFirst[l]/ratioToFirst[k]);
				
				cubicSpline_structure(s, n_dis-1, l, idx, SINGLEPOINT);

				if(ilow+1 <= n_dis-1) {
					cubicSpline_structure(s, ilow, k, idx, NORMAL);
					if(ENDTIME == OPEN_WITH_FIXED_RATIOS) s ( idx, p_index ( 0 ) );
					idx++;
				} else {
					cubicSpline_structure(s, n_dis-1, k, idx, EXTRAPOLATE);
					if(ENDTIME == OPEN_WITH_FIXED_RATIOS) s ( idx, p_index ( 0 ) );
					idx++;
				}
			}
		}
	}
	if(speedLimit) {
		for(int k=0; k<nbrCShips; k++) {
			for ( int l=0; l<1; l++ ) {
				for(int i=0; i<n_dis; i++) {
					s( idx, x_index(i, NS[k]+0 ) );
					s( idx, x_index(i, NS[k]+1 ) );
					s( idx, x_index(i, NS[k]+2 ) );
					s( idx, x_index(i, NS[k]+3 ) );
					idx++;

				}
			}
		}
	}
	return true;
}
bool OptTraj::step() {
	int iteration = folder_ptr->worhp_w.MajorIter;
	if (iteration ==0){
		if(resultfile != ""){
			ofstream eval;
			eval.open(resultfile, ios::app);
			eval << setw(12) << folder_ptr->worhp_w.NormMax_CV << " "//Constrains violation
				 << setw(12) << folder_ptr->worhp_o.F/folder_ptr->worhp_w.ScaleObj << " " ;//Constrains violation
			eval.close();
		}
	}
	std::this_thread::sleep_for (std::chrono::seconds(0));
	return true;
} 

void OptTraj::terminate() {
	int iteration = folder_ptr->worhp_w.MajorIter;
	if(resultfile != ""){
		ofstream eval;
		eval.open(resultfile, ios::app);
		eval << setw(4) << iteration <<  " "
		 << setw(12) << folder_ptr->worhp_w.NormMax_CV << " "//Constrains violation
		 << setw(12) << folder_ptr->worhp_o.F/folder_ptr->worhp_w.ScaleObj << " "; //Constrains violation
		if(folder_ptr->worhp_c.status > TerminateSuccess)
		{
			eval << setw(1) <<  "1" << endl;
			eval.close();
		}else{
			eval << setw(1) << "0" << endl;
			eval.close();
		}
	}
	if(ENDTIME == COMMON_AND_FIXED      ){
		for(int k=0; k<nbrCShips; k++) cout << "End time Ship " << k+1 << ": " << tEnd[0]*scale_t[0] << endl;
	}
	if(ENDTIME == OPEN_WITH_FIXED_RATIOS) {
		for(int k=0; k<nbrCShips; k++) cout << "End time Ship " << k+1 << ": " << ratioToFirst[k]*p(0)*scale_t[0] << endl;
	}
	if(ENDTIME == OPEN      ) {
		for(int k=0; k<nbrCShips; k++) cout << "End time Ship " << k+1 << ": " << p(k)*scale_t[k] << endl;
	}
}

void OptTraj::ToMATLAB(const std::string& filename) {

	std::ofstream of(filename);

	of.setf(std::ios::scientific);
	of << std::setprecision(9);

	of << "% TransWORHP-Result" << endl;

	of << "% ";

	for(int k=0; k<nbrCShips; k++) {
		if(ENDTIME == COMMON_AND_FIXED      ) of << std::setw(20) << tEnd[k]*scale_t[0];
		if(ENDTIME == OPEN_WITH_FIXED_RATIOS) of << std::setw(20) << p(0)*ratioToFirst[k]*scale_t[0];
		if(ENDTIME == OPEN                  ) of << std::setw(20) << p(k)*scale_t[k];
	}
	of << endl;

	for (int i = 0; i < n_dis; i++) {

		if (freetime) {
			of << std::setw(20) << (solver->T[i] * p(0) );
		} else {
			of << std::setw(20) << solver->T[i];
		}

		for(int k=0; k<nbrCShips; k++) {
			for(int j=0; j<NS[k+1]-NS[k]; j++) {
				of << std::setw(20) << x(i, j+NS[k])*scale[k][j];
			}
		}
		for(int k=0; k<nbrCShips; k++) {
			for(int j=0; j<NC[k+1]-NC[k]; j++) {
				of << std::setw(20) << u(i, j+NC[k])*scale_ctrl[k][j];
			}
		}
		of << endl;
	}
}

void OptTraj::FromMATLAB(const std::string& filename) {

	std::ifstream of(filename);
	string line;

	/// Header
	getline(of, line);

	/// Parameter p
	getline(of, line);
	string line2 = string(line,2);
	vector<double> v = ToDoubleArray(line2);

	for (int k=0; k<nbrCShips; k++) {
		if(ENDTIME == COMMON_AND_FIXED      ) /*nothing*/ ;
		if(ENDTIME == OPEN_WITH_FIXED_RATIOS) solver->X[p_index(0)] = v[0]/scale_t[0];
		if(ENDTIME == OPEN                  ) solver->X[p_index(k)] = v[k]/scale_t[k];
	}
	

	/// Data
	getline(of, line);

	vector<double> v_last;

	v = ToDoubleArray(line);
	if (freetime) v[0] = v[0]/p(0);

	v_last = v;

	for (int i=0; i<n_dis; i++) {

		while (solver->T[i] > v[0]) {

			getline(of, line);

			if (!of) {
				v_last = v;
				break;
			} else {
				v_last = v;

				v = ToDoubleArray(line);
				if (freetime) v[0] = v[0]/p(0);
			}
		}

		double vv=1;

		if (v[0]!=v_last[0])
			vv = (solver->T[i] - v_last[0]) / (v[0] - v_last[0]);

// 		for (int j=0; j<n_ode; j++) {
		for(int k=0; k<nbrCShips; k++) {
			for(int j=0; j<NS[k+1]-NS[k]; j++) {
				solver->X[x_index(i, j+NS[k])] = (v_last[j+NS[k]+1] + (v[j+NS[k]+1]-v_last[j+NS[k]+1]) * vv)/scale[k][j];

				if (i>0 && solver->twdiscretization->type==TWdiscretizationType::HermiteSimpson) {
					solver->X[x_index__(2*i-1,j+NS[k])] = (solver->X[x_index__(2*(i-1),j+NS[k])] + solver->X[x_index__(2*i,j+NS[k])])/2/scale[k][j];
				}
			}
		}

// 		for (int j=0; j<n_ctrl; j++) {
		for(int k=0; k<nbrCShips; k++) {
			for(int j=0; j<NC[k+1]-NC[k]; j++) {
				solver->X[u_index(i,j+NC[k])] = (v_last[j+NC[k]+NS[nbrCShips]+1] + (v[j+NC[k]+NS[nbrCShips]+1]-v_last[j+NC[k]+NS[nbrCShips]+1]) * vv)/scale_ctrl[k][j];

				if (i>0 && solver->twdiscretization->type==TWdiscretizationType::HermiteSimpson) {
					solver->X[u_index__(2*i-1,j+NC[k])] = (solver->X[u_index__(2*(i-1),j+NC[k]+NS[nbrCShips])] + solver->X[u_index__(2*i,j+NC[k]+NS[nbrCShips])])/2/scale_ctrl[k][j];
				}
			}
		}


	}

#ifdef TRANSWORHP_GRAPHICS
	if (solver->viewer) solver->viewer->autoScale();
#endif
}




