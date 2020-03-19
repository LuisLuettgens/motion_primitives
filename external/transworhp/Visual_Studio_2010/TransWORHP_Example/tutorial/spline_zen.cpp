/*----------------------------------------------------------------
 *
 *  Tutorial: Splineproblem
 *
 *----------------------------------------------------------------*/

#ifdef WIN32
#include "windows.h"
#endif

#include "TransWORHP.h"

using namespace std;

int MODE = 0; // ungestoert
	
class SplineZen : public TransWorhp {
public:

//	TransWorhp(s, dis, ode, ctrl, param, rand, neben, int integral=0, int zen=0);
	SplineZen(int dis) : TransWorhp("Spline",dis,3,1,0,1,0,0,3) {}
	// Zen: Letzte Zahl (3): Festlegen der Anzahl der nichtlinearen Störparameter
	// Zen: n_rand=2: Randbedingungen als Störparameter

	void zen_init(double *zen) {

		// Zen: Referenzwerte der nichtlin. Störungen
		zen[0] = (MODE==1)?.01:0;
		zen[1] = (MODE==2)?.01:0;
		zen[2] = (MODE==3)?.01:0;
	}
	
	double obj() {

		return x(n_dis-1,2);
	}

	bool obj_structure(DiffStructure &s) {

		s(0, x_index(n_dis-1,2));
		return true;
	}

	bool obj_diff(DiffStructure &s) {

		s(0, x_index(n_dis-1,2)) = 1;
		return true;
	}

	void ode(double *dx, double t, const double *x, const double *u,
			 const double *p) {

		dx[0] = x[1] + ZEN[2];
		dx[1] = u[0];
		dx[2] = u[0]*u[0];
	}

	bool ode_structure(DiffStructure &s) {

		s(0, x_indexode(1)); // dx[0] / dx[1]
		s(1, u_indexode(0)); // dx[1] / du[0]
		s(2, u_indexode(0)); // dx[2] / du[0]
		return true;
	}

	bool ode_diff(DiffStructure &s, double t, const double *x,
			      const double *u, const double *p) {

		s(0, x_indexode(1))= 1;
		s(1, u_indexode(0))= 1;
		s(2, u_indexode(0))= 2*u[0];
		return true;
	}

	void u_boundary(double *u_low, double *u_upp) {

		u_low[0] = -6;
		u_upp[0] = +6;
	}

	void x_boundary(double *x_low, double *x_upp) {

		//x_low[2] = 0;
	}

	void var_boundary(double *x_low, double *x_upp) {

		// Zen: Achtung, Box-Beschränkungen mit oberer=unterer Schranke sind (derzeit) nicht erlaubt!
		//double eps0 = .00;
		
	//	x_low[x_index(0,0)] = x_upp[x_index(0,0)] = 0 + ZEN[0];
		x_low[x_index(0,1)] = x_upp[x_index(0,1)] = 1 + ZEN[1];
		x_low[x_index(0,2)] = x_upp[x_index(0,2)] = 0;

		x_low[x_index(n_dis-1,0)] = x_upp[x_index(n_dis-1,0)] = 0;
		x_low[x_index(n_dis-1,1)] = x_upp[x_index(n_dis-1,1)] = 1;
		
		//x_low[x_index(0,0)] = 0;
	//	x_low[x_index(0,1)] = 1;
	//	x_low[x_index(0,2)] = 0;

		//x_low[x_index(n_dis-1,0)] =  0;
	  //  x_low[x_index(n_dis-1,1)] =  1;
		
		//x_upp[x_index(0,0)] = 0+ eps0 ;
	//	x_upp[x_index(0,1)] = 1 ;
	//	x_upp[x_index(0,2)] = 0 ;

		//x_upp[x_index(n_dis-1,0)] =  0 ;
		//x_upp[x_index(n_dis-1,1)] =  1 ;
		
	}

	// Zen: Beschränkungen mit Störparametern hier angeben (und nur hier, nicht in var_boundary)
	void rand(double *r) {
		r[0] = x(0,0) - ZEN[0];
	//	r[1] = x(n_dis-1,1) - ZEN[1] - 1;  // also x_{ende,1} = 1
	}
	void rand_boundary(double *r_low, double *r_upp) {
	
		for (int i=0;i<n_rand;i++) {
			r_low[i] = -.1;
			r_upp[i] = .1;
		}
	
	}
	bool rand_structure(DiffStructure &s) {
		
		s(0, x_index(0,0));		
	//	s(1, x_index(n_dis-1,1));		
		return true;
	}
	bool rand_diff(DiffStructure &s) {
		
		s(0, x_index(0,0)) = 1;	
	//	s(1, x_index(n_dis-1,1)) = 1;		
		return true;
	}
	
	void u_init(double *u, int i, int dis) {

		u[0] = -6 + (12.*i)/dis;
	}

};

///////////////////////////////////////////////////////////////////





int main(int argv, char* argc[]) {

	TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv,argc);
	TWfolder folder(&twparameter,0);

	for (int i=1;i<argv;i++) {
	
		if (strcmp(argc[i],"-m")==0) {
			i++;
			MODE = atoi(argc[i]);
		}
	}


	SplineZen ph(twparameter.NDIS);
	folder.Add(&ph);

	Viewer *viewer = 0;
	if (twparameter.PLOT) viewer = new Viewer(&twparameter);

	folder.Init();
	folder.Init(viewer);
	folder.Loop();	
	
	
	cout << "Writing matrix files..." << endl;
	TransWorhp::DoubleToMATLAB(folder.worhp_o.X, folder.worhp_o.n, "X_" + ToString(MODE) + ".m");
	TransWorhp::DoubleToMATLAB(folder.worhp_o.G, folder.worhp_o.m, "G_" + ToString(MODE) + ".m");
	TransWorhp::MatrixToMATLAB(folder.worhp_w.DF, "DF_" + ToString(MODE) + ".m");
	TransWorhp::MatrixToMATLAB(folder.worhp_w.DG, "DG_" + ToString(MODE) + ".m");
	TransWorhp::MatrixToMATLAB(folder.worhp_w.HM, "HM_" + ToString(MODE) + ".m");
	
	
 if (MODE==0) {  
	 #ifdef USE_ZEN

	 if (folder.worhp_c.ZenInit) {
	 
	char var[] = {'X','G','F','L','M',0};
	char pert[] = {'P','R','Q','B',0};
	int i=0;
	while (var[i]) {
		//os << "DX/DP" << endl;
		ofstream os( (string("ZenD") + var[i] + ".dat").c_str() );
		
		int j=0;
		while (pert[j]) {
			//os << "DX/DP" << endl;
			folder.WriteSensitivity(os, var[i], pert[j]);
			j++;
		}
		
		i++;
	}
	
	
	 double d[1000];
	 folder.GetSensitivity('X','P',1,d);
	 
	 cout << "DX1 / DP1 " << endl;
	 for (int i=0;i<twparameter.NDIS;i++) {
		 int k =ph.x_index(i,0);
		 cout << setw(12) << d[k];
	 }
	 cout <<endl;
	 
	 cout << "DX2 / DP1 " << endl;
	 for (int i=0;i<twparameter.NDIS;i++) {
		 int k = ph.x_index(i,1);
		 cout << setw(12) << d[k];
	 }
	 cout <<endl;
	 
	  cout << "DX3 / DP1 " << endl;
	 for (int i=0;i<twparameter.NDIS;i++) {
		 int k = ph.x_index(i,2);
		 cout << setw(12) << d[k];
	 }
	 cout <<endl;
	 
	  cout << "DU1 / DP1 " << endl;
	 for (int i=0;i<twparameter.NDIS;i++) {
		 int k = ph.u_index(i,0);
		 cout << setw(12) << d[k];
	 }
	 cout <<endl;
	 
	 /*new scope to define new variables */
    /*double Xnew[2], maxDP[2], maxDR[2], maxDQ[2], maxDB[2];
    double dp[2] = {0.1, 0.1}, dq[2] = {-0.2, -0.2};
    
    ZenGetMaxPert(&folder.worhp_o, &folder.worhp_w, &folder.worhp_p, &folder.worhp_c, maxDP, maxDR, maxDQ, maxDB);
    printf("\n Maximum of perturbations: \n");
    printf(" MaxDeltaP:  % .10e % .10e\n", maxDP[0], maxDP[1]);
    printf(" MaxDeltaP:  % .10e % .10e\n", maxDR[0], maxDR[1]);
    printf(" MaxDeltaP:  % .10e % .10e\n", maxDQ[0], maxDQ[1]);
    printf(" MaxDeltaP:  % .10e % .10e\n", maxDB[0], maxDB[1]);

	/*if (MODE==0) {
		TransWorhp::MatrixToMATLAB(folder.worhp_w.ZenDX, "ZenDX_" + ToString(MODE) + ".m");
		TransWorhp::MatrixToMATLAB(folder.worhp_w.ZenDMu, "ZenDMu_" + ToString(MODE) + ".m");
		TransWorhp::MatrixToMATLAB(folder.worhp_w.ZenDF, "ZenDF_" + ToString(MODE) + ".m");
		TransWorhp::MatrixToMATLAB(folder.worhp_w.ZenDF2, "ZenDF2_" + ToString(MODE) + ".m");
		TransWorhp::MatrixToMATLAB(folder.worhp_w.ZenDG, "ZenDG_" + ToString(MODE) + ".m");
	}*/
	 }


	#endif
	 
	}
	
/*	  PrintWorhpMatrix0(&w.DG);
      PrintWorhpMatrix0(&w.HM);

      PrintWorhpMatrix0(&w.ZenDX);

      PrintWorhpMatrix0(&w.ZenDF2);*/

	//delete viewer;
	
	return 0;
}

