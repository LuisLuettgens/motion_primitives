/*----------------------------------------------------------------
 *
 * Laufkatze: Einfache Trajektorie
 *
 *----------------------------------------------------------------*/

#include "laufkatze_base.h"

#include "TWfolder.h"
#include "TWproblem.h"
#include "TransWORHP.h"

double start[] = {0,0,0,0,0,5,0,0,0};
double ziel[] =  {8,0,0,0,0,4,0,0};
//double weight[] = {.015, 0.1, 0.015, 0.1, 0.1, 0.015, .1, .1};
double weight[] = {.015, 0.1, 0.015, 200, 500, 0.15, .1, 0.1};

class LaufkatzeMPC : public tw::TransWorhpProblem {
public:

	const double grav;

	XMLNode* scenenode;

	LaufkatzeMPC(const tw::TWdimension &TWdata, XMLNode *xmlmain) : TransWorhpProblem(TWdata), grav(9.81) {

		if (xmlmain) {
			scenenode = xmlmain->GetFirstChild("SCENE");
		}
	}

	void OpenWindows(tw::Viewer *viewer) override {

		viewer->ThreeD("Hallenansicht", scenenode, laufkatze3d);
	}


	void selectWindows(tw::Viewer *viewer) override {

		viewer->AddStateView(0,"Position X");
		viewer->AddStateView(1,"Geschwindigkeit X");
		viewer->AddStateView(2,"Verschiebung Last");
		viewer->AddStateView(3,"Geschwindigkeit Last");
		viewer->AddStateView(4,"Beschleunigung X");
		viewer->AddStateView(5,"Position Y");
		viewer->AddStateView(6,"Geschwindigkeit Y");
		viewer->AddStateView(7,"Beschleunigung Y");
		viewer->AddStateView(8,"u^2");
		
		viewer->AddControlView(0,"Ruck X");
		viewer->AddControlView(1,"Ruck Y");
	}

	void p_init(double *p) override {
		//p[0] = 12;
	}

	void x_init(double *x, int i, int dis) override {
		x[5] = 5;
	}
	void u_init(double *u, int i, int dis) override {
		u[0] = 0;
		u[1] = 0;
	}


	double obj() override {
		
		double ret = 0;
		
		for (int i=0;i<8;i++) {
			double tmp = x(n_dis-1,i)-ziel[i];
			ret += tmp*tmp* weight[i];
		}
			
		return ret + x(n_dis-1,8);
	}


	bool obj_structure(tw::DiffStructure &s) override {
	return false;
	s(0,x_index(n_dis-1,8) );
		return true;

	}

	bool obj_diff(tw::DiffStructure &s) override {
		return false;
		s(0,x_index(n_dis-1,8)) = 1 ;
		return true;
	}


	void ode(double *dx, double t, const double *x, const double *u, const double *p) override {

		dx[0] = x[1];
		dx[1] = x[4];
		dx[2] = x[3];
		dx[3] = (x[4] - (grav-x[7])*x[2]/x[5]);
		dx[4] = u[0];
		dx[5] = x[6];
		dx[6] = x[7];
		dx[7] = u[1];
		dx[8] = (u[0]*u[0] + u[1]*u[1]);

		//double E_pot = grav * sqrt(x[5]*x[5]-x[2]*x[2]);
		//double E_kin = (x[1]+x[3]) * (x[1]+x[3]) /2;
		
		//dx[8] += (E_pot + E_kin) * .5;
	
	}


	// Optional
	bool ode_structure(tw::DiffStructure &s) override {

		return false;
		
		s(0,x_indexode(1));
		s(0,p_indexode(0));

		s(1,x_indexode(4));
		s(1,p_indexode(0));

		s(2,x_indexode(3));
		s(2,p_indexode(0));

		s(3,x_indexode(4));
		s(3,x_indexode(7));
		s(3,x_indexode(2));
		s(3,x_indexode(5));
		s(3,p_indexode(0));

		s(4,u_indexode(0));
		s(4,p_indexode(0));

		s(5,x_indexode(6));
		s(5,p_indexode(0));

		s(6,x_indexode(7));
		s(6,p_indexode(0));

		s(7,u_indexode(1));
		s(7,p_indexode(0));

		s(8,u_indexode(0));
		s(8,u_indexode(1));
		s(8,p_indexode(0));
		
		return true;
	}

	bool ode_diff(tw::DiffStructure& s, double t, const double* x, const double* u, const double* p) override {
return false;

//return false;
		// dx[0] = x[1] * p[0];
		s(0,x_indexode(1)) = p[0];
		//s.set(0,p_indexode(0)) = x[1];


		// dx[1] = x[4] * p[0];
		s(1,x_indexode(4)) = p[0];
		//s.use(1,p_indexode(0));

		// dx[2] = x[3] * p[0];
		s(2,x_indexode(3)) = p[0];

		// dx[3] = (x[4] - (g-x[7])*x[2]/x[5])  * p[0];
		s(3,x_indexode(4)) = p[0];
		s(3,x_indexode(7)) = x[2]/x[5]* p[0];
		s(3,x_indexode(2)) = -(grav-x[7])/x[5]* p[0];
		s(3,x_indexode(5)) = (grav-x[7])*x[2]/x[5]/x[5]* p[0];

		// dx[4] = u[0] * p[0];
		s(4,u_indexode(0)) = p[0];

		// dx[5] = x[6] * p[0];
		s(5,x_indexode(6)) = p[0];

		// dx[6] = x[7] * p[0];
		s(6,x_indexode(7)) = p[0];

		// dx[7] = u[1] * p[0];
		s(7,u_indexode(1)) = p[0];

		// dx[8] = (u[0]*u[0] + u[1]*u[1])  * p[0];
		s(8,u_indexode(0)) = 2*u[0]*p[0];
		s(8,u_indexode(1)) = 2*u[1]*p[0];

		return true;
	}

	bool ode_diff_p(tw::DiffStructure& s, double t, const double *x, const double *u, const double *p, int index) override {
return false;

	
		s(0,p_indexode(0)) = x[1];
		s(1,p_indexode(0)) = x[4];
		s(2,p_indexode(0)) = x[3];
		s(3,p_indexode(0)) = x[4] - (grav-x[7])*x[2]/x[5];
		s(4,p_indexode(0)) = u[0];
		s(5,p_indexode(0)) = x[6];
		s(6,p_indexode(0)) = x[7];
		s(7,p_indexode(0)) = u[1];
		s(8,p_indexode(0)) = u[0]*u[0] + u[1]*u[1];
		
		return true;
		
	}


	void u_boundary(double *u_low, double *u_upp) override {

		u_low[0] = -1;
		u_upp[0] =  1;
		u_low[1] = -1;
		u_upp[1] =  1;

	}

	void x_boundary(double *x_low, double *x_upp) override {

		x_low[0] = 0;
		x_upp[0] = 100;

		x_low[1] = -3;
		x_upp[1] = 3;

		x_low[2] = -2;
		x_upp[2] = 2;

		x_low[3] = -10;
		x_upp[3] = 10;

		x_low[4] = -4;
		x_upp[4] = 4;

		x_low[5] = 3;//0.5;
		x_upp[5] = 6;//15;

		x_low[6] = -3;
		x_upp[6] = 3;

		x_low[7] = -10;
		x_upp[7] = 10;

		x_low[8] = 0;

	}

	void p_boundary(double *p_low, double *p_upp) override {
	
	}

	void var_boundary(double *x_low, double *x_upp) override {
	
		for (int i=0;i<9;i++) {
			x_low[x_index(0,i)] = start[i];
			x_upp[x_index(0,i)] = start[i];
		}
		
		/*for (int i=0;i<8;i++) {
			x_low[x_index(n_dis-1,i)] = ziel[i];
			x_upp[x_index(n_dis-1,i)] = ziel[i];
		}*/
		
	//	x_low[x_index(n_dis-1,1)] = x_upp[x_index(n_dis-1,1)] = 0;
		x_low[x_index(n_dis-1,2)] = x_upp[x_index(n_dis-1,2)] = 0;
		x_low[x_index(n_dis-1,3)] = x_upp[x_index(n_dis-1,3)] = 0;
		x_low[x_index(n_dis-1,4)] = x_upp[x_index(n_dis-1,4)] = 0;
	//	x_low[x_index(n_dis-1,6)] = x_upp[x_index(n_dis-1,6)] = 0;
		x_low[x_index(n_dis-1,7)] = x_upp[x_index(n_dis-1,7)] = 0;

		//x_low[x_index(0,8)] = 0;
		//x_upp[x_index(0,8)] = 0;

	}

};


/////////////////////////////////////////////////////////////////////////////



int main(int argv, char* argc[]) {

	tw::TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv,argc);
	

	std::unique_ptr<XMLNode> xml_laufkatze = tw::TWparameter::ReadParams("laufkatze.xml");


	tw::Viewer *viewer = nullptr;
	if (twparameter.PLOT) viewer = new tw::Viewer(&twparameter);

	std::ofstream mpcfile("mpc.dat");

	tw::TWdimension TWdim;
	TWdim.ID = "Laufkatze MPC";
	TWdim.n_ode = 9;
	TWdim.n_ctrl = 2;
	TWdim.n_dis = twparameter.NDIS;

	for (int K = 0; K < 60; K++) {

		tw::TWfolder folder(&twparameter, 0);

		LaufkatzeMPC ph(TWdim, xml_laufkatze.get());
		ph.setSolver(&twparameter);
		folder.Add(&ph);
		
		folder.Init();
		folder.Init(viewer);
		
		for (int i = 0 ; i < ph.n_dis; i++) {
			ph.solver->T[i] = i/(ph.n_dis-1.)*4;
		}
		
		folder.Loop(1,1);

		const int index = ph.n_dis/4-1;
		
		for (int i = 0; i < 9; i++) {
			start[i] = ph.x(index, i);
			mpcfile << start[i] << " ";
		}
		for (int i = 0; i < 2; i++) {
			mpcfile << ph.u(index,i) << " ";
		}
		mpcfile << std::endl;
		
		if(viewer) viewer->closeAll();
	}

	delete viewer;

	return 0;
}
