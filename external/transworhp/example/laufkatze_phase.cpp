/*----------------------------------------------------------------
 *
 *  Laufkatze: Bilevel-Problem
 *
 *----------------------------------------------------------------*/

#include "laufkatze_base.h"

#include "TWfolder.h"
#include "TWproblem.h"
#include "TransWORHP.h"


class LaufkatzePhase : public tw::TransWorhpProblem {
public:
	const double grav;
	XMLNode* scenenode;
	const int mode;

	LaufkatzePhase(XMLNode *xmlmain, const tw::TWdimension &TWdim, int mode) : TransWorhpProblem(TWdim),
		grav(9.81), mode(mode)
	{
		freetime = true;

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
		p[0] = 12.0;
	}

	void x_init(double *x, int i, int dis) override {
		x[5] = 5.0;
	}
	void u_init(double *u, int i, int dis) override {
		u[0] = 0.0;
		u[1] = 0.0;
	}


	double obj() override {
		return x(n_dis-1,8) + .1 * p(0);
	}


	bool obj_structure(tw::DiffStructure &s) override {
		s(0,x_index(n_dis-1,8) );
		s(0,p_index(0) );
		return true;

	}

	bool obj_diff(tw::DiffStructure &s) override {
		s(0,x_index(n_dis-1,8)) = 1.0;
		s(0,p_index(0)) = .1;
		return true;
	}


	void ode(double *dx, double t, const double *x, const double *u, const double *p) override {

		dx[0] = x[1] * p[0];
		dx[1] = x[4] * p[0];
		dx[2] = x[3] * p[0];
		dx[3] = (x[4] - (grav-x[7])*x[2]/x[5])  * p[0];
		dx[4] = u[0] * p[0];
		dx[5] = x[6] * p[0];
		dx[6] = x[7] * p[0];
		dx[7] = u[1] * p[0];
		dx[8] = (u[0]*u[0] + u[1]*u[1])  * p[0];

	}


	// Optional
	bool ode_structure(tw::DiffStructure &s) override {

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
		/*
		  for (int i=0;i<9;i++) {
		   s.use(i,p_indexode(0));
		  }
		*/
//  s.use(1,p_indexode(0)); // p_index-Variante
		//return true;
		return true;

	}

	bool ode_diff(tw::DiffStructure& s, double t, const double* x, const double* u, const double* p) override {
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
	//return false;
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

		u_low[0] = -1.0;
		u_upp[0] =  1.0;
		u_low[1] = -1.0;
		u_upp[1] =  1.0;

	}

	void x_boundary(double *x_low, double *x_upp) override {

		x_low[0] = 0.0;
		x_upp[0] = 100.0;

		x_low[1] = -3.0;
		x_upp[1] = 3.0;

		x_low[2] = -2.0;
		x_upp[2] = 2.0;

		x_low[3] = -10.0;
		x_upp[3] = 10.0;

		x_low[4] = -4.0;
		x_upp[4] = 4.0;

		x_low[5] = 0.5;
		x_upp[5] = 15.0;

		x_low[6] = -3.0;
		x_upp[6] = 3.0;

		x_low[7] = -10.0;
		x_upp[7] = 10.0;

		x_low[8] = 0.0;

	}

	void p_boundary(double *p_low, double *p_upp) override {
		
		if (mode==0) {
			p_low[0] = 1.0;
			p_upp[0] = 25.0;
		}
		else if (mode==1) {
			p_low[0] = 5.0;
			p_upp[0] = 5.0;
		}
	}

	void var_boundary(double *x_low, double *x_upp) override {

		double start[] = {0.0,0.0,0.0,0.0,0.0,5.0,0.0,0.0};
		double ziel[] =  {8.0,0.0,0.0,0.0,0.0,4.0,0.0,0.0};

		if (mode==0) {
			for (int i=0;i<8;i++) {
				x_low[x_index(0,i)] = start[i];
				x_upp[x_index(0,i)] = start[i];
			}
			for (int i=0;i<8;i++) {
				x_low[x_index(n_dis-1,i)] = ziel[i];
				x_upp[x_index(n_dis-1,i)] = ziel[i];
			}
		}
		else if (mode==1) { // nur Stopp-Trajektorie
			for (int i=1;i<5;i++) {
				x_low[x_index(n_dis-1,i)] = ziel[i];
				x_upp[x_index(n_dis-1,i)] = ziel[i];
			}
			for (int i=6;i<8;i++) {
				x_low[x_index(n_dis-1,i)] = ziel[i];
				x_upp[x_index(n_dis-1,i)] = ziel[i];
			}
					
		}
		
		x_low[x_index(0,8)] = 0.0;
		x_upp[x_index(0,8)] = 0.0;

	}

};


const int STOP = 9;
const int INDEX[] = {2,5,8,11,14,17,20,23,26};

const int STATES = 8;

class LaufkatzeFolder : public tw::TWfolder {
	
public:
	LaufkatzeFolder(tw::TWparameter *p) : TWfolder(p,STATES*STOP) {}

	void g_boundary(double *g_low, double *g_upp) override {

		for (int i=0;i<n_con;i++) {
			g_low[i] = 0;
			g_upp[i] = 0;
		}
	}

	void con(double *C) override {
		
		// zunächst nur lineare Nebenbedingungen erlaubt!
		double *X = worhp_o.X;
		
		for (int j=0;j<STOP;j++) {
			// Von der Haupt-Bahn in der Mitte...
			int index1 = phases[0]->x_index(INDEX[j],0) + phases[0]->solver->Delta1;
			
			// ...verbinden mit der Hilfbahn am Anfang
			int index2 = phases[j+1]->x_index(0,0)      + phases[j+1]->solver->Delta1;
			
			for (int i=0;i<STATES;i++) {
				//cout << j*STATES + i << endl;
				C[j*STATES + i] = X[index1+i] - X[index2+i];
				
			}
		}
	}
	
	bool con_structure(tw::DiffStructure &s) override {

		for (int j=0;j<STOP;j++) {
			int index1 = phases[0]->x_index(INDEX[j],0) + phases[0]->solver->Delta1;
			int index2 = phases[j+1]->x_index(0,0)      + phases[j+1]->solver->Delta1;
	
			for (int i=0;i<STATES;i++) {
				s(j*STATES+i, index1+i);
				s(j*STATES+i, index2+i);
			}
		}
		
		return true;
	}
	
	bool con_diff(tw::DiffStructure &s, int colindex) override {

		for (int j=0;j<STOP;j++) {
			int index1 = phases[0]->x_index(INDEX[j],0) + phases[0]->solver->Delta1;
			int index2 = phases[j+1]->x_index(0,0)      + phases[j+1]->solver->Delta1;
	
			for (int i=0;i<STATES;i++) {
				s(j*STATES+i, index1+i) = +1.0;
				s(j*STATES+i, index2+i) = -1.0;
			}
		}

		return true;
	}
	
	bool step() override {
		//cout << phases[0]->p(0) << "   " << phases[1]->p(0) << "   " << endl;
		
		//sleep(.02);
		return true;
	}
};

/////////////////////////////////////////////////////////////////////////////



int main(int argv, char* argc[]) {

	tw::TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv,argc);
	LaufkatzeFolder folder(&twparameter);
	
	// We need many discrete points for stop trajectories
	if (twparameter.NDIS<INDEX[STOP-1]+1) {
		twparameter.NDIS=INDEX[STOP-1]+1;
	}
	
	// speed up GUI
	//twparameter.showDF = false;
	//twparameter.showDG = false;
	//twparameter.showHM = false;
	
	twparameter.NDIS = 31;
	
	std::unique_ptr<XMLNode> xml_laufkatze(tw::TWparameter::ReadParams("laufkatze.xml"));
	
	tw::TWdimension TWdim;
	TWdim.ID = "Laufkatze";
	TWdim.n_dis = twparameter.NDIS;
	TWdim.n_ode = 9;
	TWdim.n_ctrl = 2;
	TWdim.n_param = 1;
	TWdim.multinode = {0,2,5,8,11,14,17,20,23,26,TWdim.n_dis-1};
	
	LaufkatzePhase ph(xml_laufkatze.get(),TWdim,0);
	//ph.LinearTimeAxis(0,1);
	
	std::vector<std::unique_ptr<LaufkatzePhase>> ph2(STOP);
	for (int i=0;i<STOP;i++) {
		tw::TWdimension twdim(TWdim);
		twdim.n_dis=11;
		ph2[i] = std::unique_ptr<LaufkatzePhase>(new LaufkatzePhase(xml_laufkatze.get(),twdim,1));
		//ph2[i]->LinearTimeAxis(0,1);
		ph2[i]->t0 = INDEX[i]/(twparameter.NDIS-1.);
	
	}
	
	ph.setSolver(&twparameter);
	folder.Add(&ph);
	for (auto & elem : ph2) {
		elem->setSolver(&twparameter);
		folder.Add(elem.get());
	}

	std::unique_ptr<tw::Viewer> viewer;
	if (twparameter.PLOT) viewer = std::unique_ptr<tw::Viewer>(new tw::Viewer(&twparameter));

	folder.Init();
	folder.Init(viewer.get());
	folder.Loop();
	
	return 0;
}
