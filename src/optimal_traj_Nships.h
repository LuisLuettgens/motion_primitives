#ifndef optimal_trag_Nships_H
#define optimal_trag_Nships_H

#ifdef WIN32
#include "windows.h"
#endif
#include <thread>       
#include <chrono> 


#include "shipviewer_base_Nships.h" 
#include "geometric_constr.h" 



class OptTraj : public tw::TransWorhpProblem {
	public:

		XMLNode* scenenode;

		vector <BaseModel*> ship_model;
		
		vector <int> NmodelS;
		vector <int> NmodelC;
		vector <int> NS;
		vector <int> NC;
		vector <vector<double>> scale;
		vector <vector<double>> scale_ctrl;
		vector<double> scale_t;
		double viewer_t_scale;
		
		int nbrCShips;
		int nbrAShips;
		int nbrAShipsMax;
		int nbrObstacles;
		int nbrObstacles_init;
		int nbrObstaclesMax;
		vector <vector<double>> var_b_init;
		vector <vector<double>> var_b_final;
		vector <vector<double>> final_s;
		
		GeometricConstr   geom_constr;
		vector<double>         distSH;
		vector<double>         distSO;
		vector<vector<double>> distSA;
		vector<vector<double>> distSS;
		vector<double> speedbounds;
		
		vector<double> tEnd;
		EndTimeBool ENDTIME;
		
		bool withCollReg;
		bool collavoidSH;
		bool collavoidSO;
		bool collavoidSA;
		bool collavoidSS;
		
		bool speedLimit;
		
		int Ninterp;
		
		vector<double> ratioToFirst;
		
		vector<double> tWgt;
		vector<double> eWgt;
		vector<vector<double>> fWgt;
		vector<vector<double>> cWgt;
		vector<double> plow;
		vector<double> pupp;
		
		vector<AISShip> AISShips;
		
		tw::TWfolder* folder_ptr;
		
		
		double *ctrl;   
		double *state; 
		double *dstate; 
		
		bool soft_final_constr;
		
		string resultfile;

		OptTraj ( XMLNode *xmlmain,
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
			);

		void OpenWindows ( tw::Viewer *viewer );
		void selectWindows ( tw::Viewer *viewer );
		void p_init ( double *p );
		double obj();
		bool obj_structure ( tw::DiffStructure &s );
		bool obj_diff ( tw::DiffStructure &s );
		void ode ( double *dx, double t, const double *x, const double *u, const double *p );
		bool ode_structure ( tw::DiffStructure &s );
		void p_boundary ( double *p_low, double *p_upp );
		void x_boundary ( double *x_low, double *x_upp );
		void u_boundary ( double *u_low, double *u_upp );
		void var_boundary ( double *x_low, double *x_upp );
		void rand(double *r);
		void rand_boundary(double *r_low, double *r_upp);
		bool rand_structure(tw::DiffStructure &s);
		bool step();
		void terminate();
		void ToMATLAB(const std::string& filename);
		void FromMATLAB(const std::string& filename);
};


#endif

