/*----------------------------------------------------------------
 *
 * Laufkatze: Regelmatrix für Störung
 *
 *----------------------------------------------------------------*/

#ifdef WIN32
#include "windows.h"
#endif

#include "TransWORHP.h"
#include "laufkatze_base.h"

using namespace std;

double start[] = {0+1,0,0,0,0,5+1,0,0};
double ziel[] =  {8,0,0,0,0,4,0,0};

#define N_STATE 8
#define N_CTRL 2

class Reference {
public:

	Reference(int ndis) : n_dis(ndis) {}

	int n_dis;
	int n_state, n_ctrl;
	double X[100][N_STATE];
	double U[100][N_CTRL];
	double T[100];

	void FromMATLAB(const std::string& filename) {

		ifstream of(filename.c_str());
		string line;

		/// Header
		getline(of, line);

		/// Parameter p
		getline(of, line);

		for (int i=0; i<n_dis; i++) {

			getline(of, line);
			vector<double> v = ToDoubleArray(line);

			T[i] = v[0];

			for (int j=0; j<N_STATE; j++) {
				X[i][j] = v[j+1];
			}

			for (int j=0; j<N_CTRL; j++) {
				U[i][j] = v[j+N_STATE+1];
			}

		}
	}
};

Reference *REF=0;
double laufkatzeplot (int &len, int &ndgl, int &nsteuer, double *t,double *x, double *u, int &i, int &index);



class LaufkatzeRegler : public TransWorhp {
public:
	double grav;
	XMLNode* scenenode;

	LaufkatzeRegler(XMLNode *xmlmain, int dis) : TransWorhp("Laufkatze", dis,9,0,18,0,0) {
	
		grav = 9.81;

		if (xmlmain)
			scenenode = xmlmain->GetFirstChild("SCENE");

	}

	void OpenWindows(Viewer *viewer) {

		viewer->ThreeD("Hallenansicht", scenenode, laufkatze3d);
		viewer->Data("u1", laufkatzeplot, 0);
		viewer->Data("u2", laufkatzeplot, 1);

	}
	
	void GetXTitle(int d, char *s) {
		if (d==0) strcpy(s, "Position X");
		if (d==1) strcpy(s, "Geschwindigkeit X");
		if (d==2) strcpy(s, "Verschiebung Last");
		if (d==3) strcpy(s, "Geschwindigkeit Last");
		if (d==4) strcpy(s, "Beschleunigung X");

		if (d==5) strcpy(s, "Position Y");
		if (d==6) strcpy(s, "Geschwindigkeit Y");
		if (d==7) strcpy(s, "Beschleunigung Y");

		if (d==8) strcpy(s, "u^2");

	}

	void p_init(double *p) {

		for (int i=0;i<n_param;i++) {
			p[i] = 0;
		}
	}

	void x_init(double *x, int i, int dis) {
		//x[5] = 5;
	}
	

	double obj() {

		double ret = 0;

		for (int i=0;i<N_STATE;i++) {
			double tmp = x(n_dis-1,i) - ziel[i];
			ret += tmp*tmp;
		}

	/*	for (int i=0;i<n_param;i++) {
			ret += p(i)*p(i);
		}*/

		ret += x(n_dis-1,8);

		return ret;
	}


	bool obj_structure(DiffStructure &s) {
		
		for (int i=0;i<N_STATE;i++) {
			s(0,x_index(n_dis-1,i) );
		}
		
		//for (int i=0;i<n_param;i++) {
		//	s(0,p_index(i) );
		//}
		
		
		s(0,x_index(n_dis-1,8) );
		
		return true;
	}

	bool obj_diff(DiffStructure &s) {
		
		for (int i=0;i<N_STATE;i++) {
			s(0,x_index(n_dis-1,i) ) = 2 * (x(n_dis-1,i) - ziel[i]);
		}

		//for (int i=0;i<n_param;i++) {
		//	s(0,p_index(i) ) = 2 * p(i);
		//}
		
		s(0,x_index(n_dis-1,8) ) = 1;
		
		return true;

	}

	/* i=0..9, j=0,1 */
	double K(int i, int j) {return p(i*N_CTRL+j);}

	double CalcU(int dis, int ctrl) {

		double uu = 0;

		for (int i=0;i<N_STATE;i++) {
			uu += K(i,ctrl) * (x(dis,i) - REF->X[dis][i]);
		}

		return uu;
	}


	void ode(double *dx, double t, const double *x, const double *u, const double *p) {

		int dis = (int) (t / T[n_dis-1] * (n_dis-1) + 0.5);

		double u0 = REF->U[dis][0] + CalcU(dis,0);
		double u1 = REF->U[dis][1] + CalcU(dis,1);

		dx[0] = x[1];
		dx[1] = x[4];
		dx[2] = x[3];
		dx[3] = (x[4] - (grav-x[7])*x[2]/x[5]);
		dx[4] = u0;
		dx[5] = x[6];
		dx[6] = x[7];
		dx[7] = u1;
		dx[8] = CalcU(dis,0)*CalcU(dis,0) + CalcU(dis,1)*CalcU(dis,1);

	}

	void x_boundary(double *x_low, double *x_upp) {

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

	void p_boundary(double *p_low, double *p_upp) {

		for (int i=0;i<n_param;i++) {

			p_low[i] = -2.0;
			p_upp[i] = 2.0;
		}
	}

	void var_boundary(double *x_low, double *x_upp) {

		for (int i=0;i<N_STATE;i++) {
			x_low[x_index(0,i)] = start[i];
			x_upp[x_index(0,i)] = start[i];
		}

		x_low[x_index(0,8)] = 0;
		x_upp[x_index(0,8)] = 0;

	}

	void terminate() {

		cout << "Matrix K:" << endl;
		cout.setf(ios::fixed);
		for (int j=0;j<N_CTRL;j++) {
			for (int i=0;i<N_STATE;i++) {
				cout << setprecision(6) << setw(14) << K(i,j);

			}
			cout << endl;
		}
		
		ofstream os("K.dat");
		for (int j=0;j<N_CTRL;j++) {
			for (int i=0;i<N_STATE;i++) {
				os << setprecision(6) << setw(14) << K(i,j);

			}
			os << endl;
		}
		

	}

};



LaufkatzeRegler *PH=0;

double laufkatzeplot (int &len, int &ndgl, int &nsteuer, double *t,double *x, double *u, int &i, int &index) {

	if (index==0) return REF->U[i][0] + PH->CalcU(i,index);
	if (index==1) return REF->U[i][1] + PH->CalcU(i,index);

}


/////////////////////////////////////////////////////////////////////////////



int main(int argv, char* argc[]) {

	TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv,argc);
	TWfolder folder(&twparameter,0);

	XMLNode *xml_laufkatze = TWparameter::ReadParams("laufkatze.xml");


	Reference ref(twparameter.NDIS);
	ref.FromMATLAB("referenz.dat");
	REF = &ref;


	LaufkatzeRegler ph(xml_laufkatze, twparameter.NDIS);
	folder.Add(&ph);
	PH = &ph;

	Viewer *viewer = 0;
	if (twparameter.PLOT) viewer = new Viewer(&twparameter);

	folder.Init();
	folder.Init(viewer);

	for (int i=0;i<ref.n_dis;i++) ph.T[i] = ref.T[i];
	

	folder.Loop(2);

	delete viewer;
	delete xml_laufkatze;

	return 0;
}

