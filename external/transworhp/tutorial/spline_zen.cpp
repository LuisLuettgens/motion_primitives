/*----------------------------------------------------------------
 *
 *  Tutorial: Splineproblem mit ZEN
 *
 *----------------------------------------------------------------*/

#include "TWfolder.h"
#include "TWproblem.h"
#include "TransWORHP.h"

#include <cstring>

int MODE = 0; // ungestoert

class SplineZen : public tw::TransWorhpProblem {
public:

	SplineZen(int dis) : TransWorhpProblem(tw::TWdimension("Spline",dis,3,1,0,1,0,0,3)) {}
	// Zen: Letzte Zahl (3): Festlegen der Anzahl der nichtlinearen Stoerparameter
	// Zen: n_rand=1: Randbedingungen als Stoerparameter

	void zen_init(double *zen) {

		// Zen: Referenzwerte der nichtlin. Stoerungen
		zen[0] = (MODE==1)? 0.01 : 0.0;
		zen[1] = (MODE==2)? 0.01 : 0.0;
		zen[2] = (MODE==3)? 0.01 : 0.0;
	}

	double obj() override {

		return x(n_dis-1,2);
	}

	bool obj_structure(tw::DiffStructure &s) override {

		s(0, x_index(n_dis-1,2));
		return true;
	}

	bool obj_diff(tw::DiffStructure &s) override {

		s(0, x_index(n_dis-1,2)) = 1;
		return true;
	}

	void ode(double *dx, double t, const double *x, const double *u, const double *p) override {

		dx[0] = x[1] + solver->ZEN[2];
		dx[1] = u[0];
		dx[2] = u[0]*u[0];
	}

	bool ode_structure(tw::DiffStructure &s) override {

		s(0, x_indexode(1)); // dx[0] / dx[1]
		s(1, u_indexode(0)); // dx[1] / du[0]
		s(2, u_indexode(0)); // dx[2] / du[0]
		return true;
	}

	bool ode_diff(tw::DiffStructure &s, double t, const double *x, const double *u, const double *p) override {

		s(0, x_indexode(1))= 1;
		s(1, u_indexode(0))= 1;
		s(2, u_indexode(0))= 2*u[0];
		return true;
	}

	void u_init(double *u, int i, int dis) override {
		u[0] = -6 + (12.*i)/dis;
	}

	void u_boundary(double *u_low, double *u_upp) override {

		u_low[0] = -6;
		u_upp[0] = +6;
	}

	void x_boundary(double *x_low, double *x_upp) override {

		x_low[2] = 0;
	}

	void var_boundary(double *x_low, double *x_upp) override {

		// Zen: Achtung, Box-Beschraenkungen mit oberer=unterer Schranke sind (derzeit) nicht erlaubt!
		//double eps0 = .00;
		
	//	x_low[x_index(0,0)] = x_upp[x_index(0,0)] = 0 + solver->ZEN[0];
		x_low[x_index(0,1)] = x_upp[x_index(0,1)] = 1 + solver->ZEN[1]; // hat keine Sens.!
		x_low[x_index(0,2)] = x_upp[x_index(0,2)] = 0;

		x_low[x_index(n_dis-1,0)] = x_upp[x_index(n_dis-1,0)] = 0;
		x_low[x_index(n_dis-1,1)] = x_upp[x_index(n_dis-1,1)] = 1;
	}

	// Zen: Beschraenkungen mit Stoerparametern hier angeben (und nur hier, nicht in var_boundary)
	// solch eine einfache Beschr. wird aber auch durch Box-B. automatisch
	// berechnet.
	void rand(double *r) override {
		r[0] = x(0,0) - solver->ZEN[0];
	//	r[1] = x(n_dis-1,1) - solver->ZEN[1] - 1;  // also x_{ende,1} = 1
	}

	void rand_boundary(double *r_low, double *r_upp) override {
	
		for (int i = 0; i < n_rand; i++) {
			r_low[i] = -.1;
			r_upp[i] = .1;
		}
	}

	bool rand_structure(tw::DiffStructure &s) override {
		
		s(0, x_index(0,0));
	//	s(1, x_index(n_dis-1,1));
		return true;
	}

	bool rand_diff(tw::DiffStructure &s) override {
		
		s(0, x_index(0,0)) = 1;
	//	s(1, x_index(n_dis-1,1)) = 1;
		return true;
	}
};

///////////////////////////////////////////////////////////////////


int main(int argv, char* argc[]) {

	tw::TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv, argc);
	tw::TWfolder folder(&twparameter, 0);

	for (int i = 1; i < argv; i++) {
		if (strcmp(argc[i], "-m") == 0) {
			i++;
			MODE = atoi(argc[i]);
			std::cout << "MODE = " << MODE << std::endl;
		}
	}

	SplineZen ph(twparameter.NDIS);
	ph.setSolver(&twparameter);
	folder.Add(&ph);

	tw::Viewer *viewer = nullptr;
	if (twparameter.PLOT) viewer = new tw::Viewer(&twparameter);

	folder.Init();
	folder.Init(viewer);
	folder.Loop();

	std::cout << "Writing matrix files..." << std::endl;
	tw::TransWorhp::DoubleToMATLAB(folder.worhp_o.X, folder.worhp_o.n, "X_" + std::to_string(MODE) + ".m");
	tw::TransWorhp::DoubleToMATLAB(folder.worhp_o.G, folder.worhp_o.m, "G_" + std::to_string(MODE) + ".m");
	tw::TransWorhp::MatrixToMATLAB(folder.worhp_w.DF, "DF_" + std::to_string(MODE) + ".m");
	tw::TransWorhp::MatrixToMATLAB(folder.worhp_w.DG, "DG_" + std::to_string(MODE) + ".m");
	tw::TransWorhp::MatrixToMATLAB(folder.worhp_w.HM, "HM_" + std::to_string(MODE) + ".m");

	if (MODE==0) {

		if (folder.worhp_c.ZenInit) {
		
		// komplette Sensitivitaeten
		const char var[] = {'X','G','F','L','M',0};
		const char pert[] = {'P','R','Q','B',0};
		int i=0;
		while (var[i]) {
			//os << "DX/DP" << endl;
			std::ofstream os( (std::string("ZenD") + var[i] + ".dat").c_str() );
			
			int j = 0;
			while (pert[j]) {
				//os << "DX/DP" << endl;
				folder.WriteSensitivity(os, var[i], pert[j]);
				j++;
			}
			
			i++;
		}
		

		// Sortierung nach Zustaenden
		double d[1000];

		// Sensitivitaet bzgl ZEN[0] (0, nicht 1, wegen FORTRAN!!!)
		folder.GetSensitivity('X', 'P', 1, d);
		
		// Sens. von erstem Zustand?
		std::cout << "DX1 / DP1 " << std::endl;
		for (int i = 0; i < twparameter.NDIS; i++) {
			const int k = ph.x_index(i,0);
			std::cout << std::setw(12) << d[k];
		}
		std::cout << std::endl;
		
		std::cout << "DX2 / DP1 " << std::endl;
		for (int i = 0; i < twparameter.NDIS; i++) {
			const int k = ph.x_index(i,1);
			std::cout << std::setw(12) << d[k];
		}
		std::cout << std::endl;
		
		std::cout << "DX3 / DP1 " << std::endl;
		
		for (int i = 0; i < twparameter.NDIS; i++) {
			const int k = ph.x_index(i,2);
			std::cout << std::setw(12) << d[k];
		}
		std::cout << std::endl;
		
		std::cout << "DU1 / DP1 " << std::endl;
		for (int i = 0; i < twparameter.NDIS; i++) {
			const int k = ph.u_index(i,0);
			std::cout << std::setw(12) << d[k];
		}
		std::cout << std::endl;
		
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
	}

	/*PrintWorhpMatrix0(&w.DG);
	PrintWorhpMatrix0(&w.HM);
	PrintWorhpMatrix0(&w.ZenDX);
	PrintWorhpMatrix0(&w.ZenDF2);*/

	delete viewer;

	return 0;
}
