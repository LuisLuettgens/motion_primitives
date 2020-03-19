/*----------------------------------------------------------------
 *
 * Hang Glider aus Betts S.283ff
 * maximiere Reichweite
 *
 *----------------------------------------------------------------*/

#include "TWfolder.h"
#include "TWproblem.h"
#include "TransWORHP.h"

class HangGliderPhase : public tw::TransWorhpProblem {
public:

	const double uM, R, C0, k, m, S, rho, g;
	HangGliderPhase(int dis) : TransWorhpProblem(tw::TWdimension("Hang Glider",dis,4,1,1,0,0)),
		uM(2.5),
		R(100.0),
		C0(0.034),
		k(0.069662),
		m(100.0),
		S(14.0),
		rho(1.13),
		g(9.80665) {
	
		freetime = true;
	}
	

	void x_init(double *x, int i, int dis) override {
		/*
		x[0] = 0;
		x[1] = 1000;
		x[2] = 13.23;
		x[3] = -1.29;
		*/
		
		x[0] = 1000-1000*(dis-i)/dis;
		x[1] = 100*(dis-i)/dis+900;
		x[2] = 13.23;
		x[3] = -1.29;
	}

	void u_init(double *u, int i, int dis) override {

		u[0] = 1.0;
	}
	
	void p_init(double *p) override {
		p[0] = 100.0;
	}

	double obj() override {
		return -x(n_dis-1, 0);
	}


	bool obj_structure(tw::DiffStructure &s) override {

		s(0,x_index(n_dis-1,0));
		return true;
	}

	bool obj_diff(tw::DiffStructure &s) override {

		s(0,x_index(n_dis-1,0)) = -1.0;
		return false;
	}
	
	void ode(double *dx, double t, const double *x, const double *u, const double *p) override {

		double CD;
		double D, L, X, ua, Vy, vr, sinn, cosn;

		CD = C0+k*u[0]*u[0];
		
		
		X = (x[0]/R-2.5)*(x[0]/R-2.5);
		
		ua = uM*(1-X)*exp(-X);
		
		Vy = x[3] - ua;
		
		vr = sqrt(x[2]*x[2] + Vy*Vy);
		
		sinn = Vy/vr;
		cosn = x[2]/vr;
		
		D = 0.5*CD*rho*S*vr*vr;
		
		L = 0.5*u[0]*rho*S*vr*vr;
		
	  
		dx[0] = x[2];
		dx[1] = x[3];
		dx[2] = 1/m*( -L * sinn -D * cosn );
		dx[3] = 1/m*( L * cosn - D * sinn -m*g );
		
		dx[0] *= p[0];
		dx[1] *= p[0];
		dx[2] *= p[0];
		dx[3] *= p[0];
		
		//cout << dx[0] << " " << dx[1] << " " << dx[2] << " " << dx[3] << endl;
	}

	bool ode_structure(tw::DiffStructure &s) override {
		
		
		return false;
	}

	bool ode_diff(tw::DiffStructure &s, double t, const double *x, const double *u, const double *p) override {
		
		
		return false;
	}
	
	bool ode_diff_p(tw::DiffStructure& s, double t, const double *x, const double *u, const double *p, int index) override {
	
		return false;
	}

	void u_boundary(double *u_low, double *u_upp) override {
		u_low[0] = 0.0;
		u_upp[0] = 1.4;
	}

	void x_boundary(double *x_low, double *x_upp) override {

		x_low[0] = 0.0;

	}
	
	void p_boundary(double *p_low, double *p_upp) override {
		p_low[0] = 80.0;
		p_upp[0] = 120.0;
	}

	void var_boundary(double *x_low, double *x_upp) override {

		x_low[x_index(0,0)] = x_upp[x_index(0,0)] = 0.0;
		x_low[x_index(0,1)] = x_upp[x_index(0,1)] = 1000.0;
		x_low[x_index(0,2)] = x_upp[x_index(0,2)]= 13.227567500;
		x_low[x_index(0,3)] = x_upp[x_index(0,3)] = -1.2875005200;
		
		//x_low[x_index(n_dis-1,0)] = x_upp[x_index(n_dis-1,0)] = ; //frei
		x_low[x_index(n_dis-1,1)] = x_upp[x_index(n_dis-1,1)] = 900.0;
		x_low[x_index(n_dis-1,2)] = x_upp[x_index(n_dis-1,2)] = 13.227567500;
		x_low[x_index(n_dis-1,3)] = x_upp[x_index(n_dis-1,3)] = -1.2875005200;
	}
};

/////////////////////////////////////////////////////////////////////////////

int main(int argv, char* argc[]) {

	tw::TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv, argc);

	tw::Viewer *viewer = nullptr;

	if (twparameter.PLOT) viewer = new tw::Viewer(&twparameter);

	tw::TWfolder folder(&twparameter, 0);
	
	HangGliderPhase ph(twparameter.NDIS);
	ph.setSolver(&twparameter);
	folder.Add(&ph);

	folder.Init();
	folder.Init(viewer);
	folder.Loop();

	delete viewer;

	return 0;
}
