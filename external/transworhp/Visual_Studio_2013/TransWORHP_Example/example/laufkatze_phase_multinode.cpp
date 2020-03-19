/*----------------------------------------------------------------
 *
 *  Laufkatze: Bilevel-Problem
 *
 *----------------------------------------------------------------*/

#ifdef WIN32
#include "windows.h"
#endif

#include "ExplTransWORHP.h"
#include "laufkatze_base.h"

using namespace std;


class LaufkatzePhase : public ExplTransWorhp {
public:
	double grav;
	XMLNode* scenenode;
	int mode;

	LaufkatzePhase(XMLNode *xmlmain, int dis, vector<int> &multinode, int mode) : ExplTransWorhp("Laufkatze", dis, multinode, 9,2,1,0,0), mode(mode) {

		freetime = 1;
		grav = 9.81;

		if (xmlmain)
			scenenode = xmlmain->GetFirstChild("SCENE");
	}
	
	void selectWindows(Viewer *viewer) override {

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
		p[0] = 12;
	}

	void x_init(double *x, int i, int dis) override {
		x[5] = 5;
	}

	void u_init(double *u, int i, int dis) override {
		u[0] = 0;
		u[1] = 0;
	}

	double obj() override {
		return x(n_dis-1,8) + .1 * p(0);
	}

	bool obj_structure(DiffStructure &s) override {
		s(0,x_index(n_dis-1,8) );
		s(0,p_index(0) );
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
	
	bool ode_diff2(DiffStructure& s, double t, const double* x, const double* u, const double* p) override {
		
		s(0,x_indexode(1)) = p[0];
		
		s(1,x_indexode(4)) = p[0];
		
		s(2,x_indexode(3)) = p[0];
		
		s(3,x_indexode(4)) = p[0];
		s(3,x_indexode(7)) = x[2]/x[5]* p[0];
		s(3,x_indexode(2)) = -(grav-x[7])/x[5]* p[0];
		s(3,x_indexode(5)) = (grav-x[7])*x[2]/x[5]/x[5]* p[0];
		
		s(5,x_indexode(6)) = p[0];
		
		s(6,x_indexode(7)) = p[0];
		
		return true;
	}

	// Optional
	bool ode_structure(DiffStructure &s) override {

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

		x_low[5] = 0.5;
		x_upp[5] = 15;

		x_low[6] = -3;
		x_upp[6] = 3;

		x_low[7] = -10;
		x_upp[7] = 10;

		x_low[8] = 0;
	}

	void p_boundary(double *p_low, double *p_upp) override {
	  
		if (mode==0) {
			p_low[0] = 1;
			p_upp[0] = 25;
		}
		else if (mode==1) {
			p_low[0] = 5;
			p_upp[0] = 5;
		}
	}

	void var_boundary(double *x_low, double *x_upp) override {

		double start[] = {0,0,0,0,0,5,0,0};
		double ziel[] =  {8,0,0,0,0,4,0,0};

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

		x_low[x_index(0,8)] = 0;
		x_upp[x_index(0,8)] = 0;
	}

	void neben_boundary(double *c_low, double *c_upp) override {

		/*c_low[0] = 0; //-65./180.*M_PI;
		c_upp[0] = 100;

		c_low[1] = -3; //-65./180.*M_PI;
		c_upp[1] = 3;

		c_low[2] = -3; //-65./180.*M_PI;
		c_upp[2] = 3;

		c_low[3] = -3; //-65./180.*M_PI;
		c_upp[3] = 3;*/
	}

	void neben(double *c, double t, const double *x, const double *u, const double *p) override {
		/*c[0] = x[0];
		c[1] = x[1];
		c[2] = x[2];
		c[3] = x[3];	*/
	}

	bool neben_structure(DiffStructure &s) override {
		return false;
		s(0,x_indexode(0));
		s(1,x_indexode(1));
		s(2,x_indexode(2));
		s(3,x_indexode(3));
		return true;
	}
};



/////////////////////////////////////////////////////////////////////////////
// Phasen

/*
* die Werte in INDEX machen bei expl-TW Probleme, da dann auf flasche Inezies zugegriffen wird!!
*/

const int STOP = 9;
int INDEX[] = {2,5,8,11,14,17,20,23,26};

const int STATES = 8;

class LaufkatzeFolder : public TWfolder {

	public:
		LaufkatzeFolder(TWparameter *p) : TWfolder(p,STATES*STOP) {}

	void g_boundary(double *x_low, double *x_upp) override {
	 
		for (int i=0;i<n_con;i++) {
			x_low[i] = 0;
			x_upp[i] = 0;
		}
	 }

	void con(double *C) override {
		
		// zunaechst nur lineare Nebenbedingungen erlaubt!
		double *X = worhp_o.X;
		
		for (int j=0;j<STOP;j++) {
			// Von der Haupt-Bahn in der Mitte...
			int index1 = phases[0]->x_index(INDEX[j],0) + phases[0]->Delta1;
			
			// ...verbinden mit der Hilfbahn am Anfang
			int index2 = phases[j+1]->x_index(0,0)                  + phases[j+1]->Delta1;
			
			for (int i=0;i<STATES;i++) {
				//cout << j*STATES + i << endl;
				C[j*STATES + i] = X[index1+i] - X[index2+i];
				
			}
		}	
	}
	
	bool con_structure(DiffStructure &s) override {

		for (int j=0;j<STOP;j++) {
			int index1 = phases[0]->x_index(INDEX[j],0) + phases[0]->Delta1;
			int index2 = phases[j+1]->x_index(0,0)                  + phases[j+1]->Delta1;
	
			for (int i=0;i<STATES;i++) {
				s(j*STATES+i, index1+i);
				s(j*STATES+i, index2+i);
			}
		}
			
		return true;
	}
	
	bool con_diff(DiffStructure &s, int colindex) override {

		for (int j=0;j<STOP;j++) {
			int index1 = phases[0]->x_index(INDEX[j],0) + phases[0]->Delta1;
			int index2 = phases[j+1]->x_index(0,0)                  + phases[j+1]->Delta1;
	
			for (int i=0;i<STATES;i++) {
				s(j*STATES+i, index1+i) = +1;
				s(j*STATES+i, index2+i) = -1;
			}
		}

		return true;
	}
	 
	bool step() override {

		return true;
	}
};
/////////////////////////////////////////////////////////////////////////////



int main(int argv, char* argc[]) {

	TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv,argc);
	
	if (twparameter.NDIS<INDEX[STOP-1]+1) {
		twparameter.NDIS=INDEX[STOP-1]+1;
	}
	
	twparameter.NDIS = 31;

	LaufkatzeFolder folder(&twparameter);


	twparameter.twdiscretization = TWdiscretization(TW_MultipleShooting,0,0);
	twparameter.USERHM=-1;
	twparameter.USERDF=-1;
	//twparameter.USERDG=-1;
	//twparameter.hessianstructure = 1; //!!!!!!!


	XMLNode *xml_laufkatze = TWparameter::ReadParams("laufkatze.xml");

	vector<int> multinode1, multinode2, multinode3;
	multinode1.push_back(0);
	multinode1.push_back(twparameter.NDIS*.25);
	multinode1.push_back(twparameter.NDIS*.75);
	multinode1.push_back(twparameter.NDIS-1);

	multinode2.push_back(0);
	multinode2.push_back(twparameter.NDIS*.5);
	multinode2.push_back(twparameter.NDIS-1);
	
	int INDEX[] = {2,5,8,11,14,17,20,23,26};
	
	multinode3.push_back(0);
	for (auto &elem : INDEX) {
		multinode3.push_back(elem);
	}
	multinode3.push_back(twparameter.NDIS-1);
	
	
	//LaufkatzePhase ph(xml_laufkatze, twparameter.NDIS, multinode1);
	//LaufkatzePhase ph2(xml_laufkatze,twparameter.NDIS*2, multinode2);
	//ph2.t0 = 2/(twparameter.NDIS-1.);
	
	LaufkatzePhase ph(xml_laufkatze, twparameter.NDIS, multinode3,0);
	
	LaufkatzePhase ph2(xml_laufkatze,twparameter.NDIS, multinode3,1);
	LaufkatzePhase ph3(xml_laufkatze,twparameter.NDIS, multinode3,1);
	LaufkatzePhase ph4(xml_laufkatze,twparameter.NDIS, multinode3,1);
	LaufkatzePhase ph5(xml_laufkatze,twparameter.NDIS, multinode3,1);
	LaufkatzePhase ph6(xml_laufkatze,twparameter.NDIS, multinode3,1);
	LaufkatzePhase ph7(xml_laufkatze,twparameter.NDIS, multinode3,1);
	LaufkatzePhase ph8(xml_laufkatze,twparameter.NDIS, multinode3,1);
	LaufkatzePhase ph9(xml_laufkatze,twparameter.NDIS, multinode3,1);
	LaufkatzePhase ph10(xml_laufkatze,twparameter.NDIS, multinode3,1);
	
	ph2.t0 = INDEX[0]/(twparameter.NDIS-1.);
	ph3.t0 = INDEX[1]/(twparameter.NDIS-1.);
	ph4.t0 = INDEX[2]/(twparameter.NDIS-1.);
	ph5.t0 = INDEX[3]/(twparameter.NDIS-1.);
	ph6.t0 = INDEX[4]/(twparameter.NDIS-1.);
	ph7.t0 = INDEX[5]/(twparameter.NDIS-1.);
	ph8.t0 = INDEX[6]/(twparameter.NDIS-1.);
	ph9.t0 = INDEX[7]/(twparameter.NDIS-1.);
	ph10.t0 = INDEX[8]/(twparameter.NDIS-1.);
	
	ph.butcherInit(twparameter.butchertableau,twparameter.stepsize );
	ph2.butcherInit(twparameter.butchertableau,twparameter.stepsize );
	ph3.butcherInit(twparameter.butchertableau,twparameter.stepsize );
	ph4.butcherInit(twparameter.butchertableau,twparameter.stepsize );
	ph5.butcherInit(twparameter.butchertableau,twparameter.stepsize );
	ph6.butcherInit(twparameter.butchertableau,twparameter.stepsize );
	ph7.butcherInit(twparameter.butchertableau,twparameter.stepsize );
	ph8.butcherInit(twparameter.butchertableau,twparameter.stepsize );
	ph9.butcherInit(twparameter.butchertableau,twparameter.stepsize );
	ph10.butcherInit(twparameter.butchertableau,twparameter.stepsize );
	
	folder.Add(&ph);
	folder.Add(&ph2);
	folder.Add(&ph3);
	folder.Add(&ph4);
	folder.Add(&ph5);
	folder.Add(&ph6);
	folder.Add(&ph7);
	folder.Add(&ph8);
	folder.Add(&ph9);
	folder.Add(&ph10);


	Viewer *viewer = nullptr;
	if (twparameter.PLOT) viewer = new Viewer(&twparameter);

	folder.Init();
	folder.Init(viewer);

	folder.Loop();

	delete viewer;
	delete xml_laufkatze;

	return 0;
}

