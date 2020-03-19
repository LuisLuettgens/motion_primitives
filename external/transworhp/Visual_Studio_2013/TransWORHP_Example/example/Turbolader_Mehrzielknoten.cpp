/*----------------------------------------------------------------
 *
 * Parameteridentifikation
 *
 *----------------------------------------------------------------*/

#ifdef WIN32
#include "windows.h"
#endif

#define _USE_MATH_DEFINES
#include "math.h"
#include "ExplTransWORHP.h"
#include "TransWORHP.h"


// Hier kann eingestellt werden ob Explizit oder mit voller Diskretisierung gerechnet werden soll
#define _USE_EXPL_TransWORHP
//#define _USE_TransWORHP

using namespace std;

double J = 5.0E-5;
double cp_FL = 1012.0;
double cp_AL = 1063.0;
double k_FL = 1.4;
double k_AL = 1.37;
double skal = 3750.0;
double dT = 0.01;

double *test2 = new double[10];
double *Parameter = new double[10];


#ifdef _USE_EXPL_TransWORHP
class Turbo_Expl : public ExplTransWorhp {
  #endif
  #ifdef _USE_TransWORHP
class Turbo_Expl : public TransWorhp {
  #endif
public:
	// Structure in dem alle Messwerte als Arrays gespeichert werden
	struct constMeasuredData {
		double* m_vrd;
		double* T_vVrd;
		double* p_nVrd;
		double* p_vVrd;
		double* m_Trb;
		double* T_vTrb;
		double* p_nTrb;
		double* p_vTrb;
		double* sVtg;
		double* nATL;
	};


	// Anlegen spezieller Structures fuer die Messdaten und Parameter (Kennlinien und Kennfelder)
	struct constMeasuredData Data;
	// time Vektor, speziell fuer die Iterpolation
	double *time;
	


#ifdef _USE_EXPL_TransWORHP	
	Turbo_Expl( int dis, vector<int> &multinode) : ExplTransWorhp("Turbo_Expl",dis,multinode,1,0,6,0,0) {
#endif
#ifdef _USE_TransWORHP	
	Turbo_Expl( int dis, vector<int> &multinode) : TransWorhp("Turbo_Full",dis,1,0,6,0,0) {
#endif
		// Messdaten einlesen
		double x,y,z,p,q,r,s,t,u; 
		//Dynamischer Speicher in der Groesse der einzulesenden Messdaten wird angelegt
		Data.m_vrd = new double[dis];
		Data.T_vVrd = new double[dis];
		Data.p_nVrd = new double[dis];
		Data.p_vVrd = new double[dis];
		Data.m_Trb = new double[dis];
		Data.T_vTrb = new double[dis];
		Data.p_nTrb = new double[dis];
		Data.p_vTrb = new double[dis];
		Data.sVtg = new double[dis];
		Data.nATL = new double[dis];
		time = new double[dis];
		
		
		//oeffnen der Messdatei
		ifstream in("Data_In_Turbolader.txt",ifstream::in);
		if (!in.good()){
			std::cout<<endl<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Fehler beim lesen der Messdatei!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
			exit(1);
		}else{
			std::cout<<"Eingangsdaten eingelesen!"<<endl;
		}
		int i=0;
		// Schleife die Zeile fuer Zeile der Messdatei einliest
		while(in){	
			
			in >> x >> y >> z >> p >> q >> r >> s >> t >> u ;
			if(!in) break;
			if(i<dis){
			  Data.m_vrd[i] = x;
			  Data.T_vVrd[i] = y;
			  Data.p_nVrd[i] = z;
			  Data.p_vVrd[i] = p;
			  Data.m_Trb[i] = q;
			  Data.T_vTrb[i] = r;
			  Data.p_nTrb[i] = s;
			  Data.p_vTrb[i] = t;
			  Data.sVtg[i] = u;
			  time[i] = i*dT;
			  i++;
			}
		}
		in.close();
		//oeffnen der Messdatei
		ifstream in2("Data_Out_Turbolader.txt",ifstream::in);
		if (!in2.good()){
			std::cout<<endl<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Fehler beim lesen der Messdatei!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
			exit(1);
		}else{
			std::cout<<"Ausgangssdaten eingelesen!"<<endl;
		}
		i=0;
		// Schleife die Zeile fuer Zeile der Messdatei einliest
		while(in2){	
			
			in2 >> x ;
			if(!in2) break;
			if(i<dis){
			  Data.nATL[i] = x;
			  i++;
			}
		}
		in2.close();
		
	}
	
	void OpenWindows(Viewer *viewer) override {
		int j=1;
		// Viewer wird konstruiert, 4 Zustaendes deshalb 4 Plots, in jeden Plot wird eine CompareCurfe aus Messdaten gezeichnet
		//viewer->plots[0]->AddCompareCurve2(T,Data.nATL_skal,j,  n_dis);
		viewer->plots[0]->AddCompareCurve2(T,Data.nATL,j,  n_dis);
	}
	
	

	void p_init(double *p) override {
	  
	  #ifdef _USE_TransWORHP
		p[0] = 0.1;
		p[1] = 0.1;
		p[2] = 0.1;
		p[3] = 0.1;
		p[4] = 0.1;
		p[5] = 0.1;
	  #endif
	  #ifdef _USE_EXPL_TransWORHP
		p[0] = 0.0;
		p[1] = 0.0;
		p[2] = 0.0;
		p[3] = 0.0;
		p[4] = 0.0;
		p[5] = 0.0;
	  #endif
		
	}

	void x_init(double *x, int i, int dis) override {
		x[0] = 0.4259;
	}



	double obj() override {
	  double sum = 0;
	  for(int i=0;i<n_dis;i++){
	      sum += pow(x(i,0) - Data.nATL[i],2);
	  }
	  return sum/n_dis;
	}

	#ifdef _USE_TransWORHP
	bool obj_structure(DiffStructure &s) {
	  
	  	
	  for(int i=0;i<n_dis;i++){
	    s(0,x_index(i,0) );
	  }
		
		return true;


	}

	bool obj_diff(DiffStructure &s) {
	
		for(int i=0;i<n_dis;i++){
		  s(0,x_index(i,0) ) = 2.0*(x(i,0) - Data.nATL[i])/n_dis;
		 }
		
		return true;
	   
	}
	#endif


	void ode(double *dx, double t, const double *x, const double *u, const double *p) override {
	  
	double w_t,w_v,P_V,P_T,m_V,T_vV,p_nV, p_vV,m_T,T_vT,p_nT,p_vT,sVtg;  
	  
	double x_skal = x[0]*skal;
	
	//Messdaten auswerten
	getParams(t,Parameter);
	
	m_V = Parameter[0];
	T_vV = Parameter[1];
	p_nV = Parameter[2];
	p_vV = Parameter[3];
	m_T = Parameter[4];
	T_vT = Parameter[5];
	p_nT = Parameter[6];
	p_vT = Parameter[7];
	sVtg = Parameter[8];
		
	
	//Auswertung der Polynome
	w_v = p[0] + p[1]*x[0] + p[2]*p_vV/p_nV;
	w_t = p[3] + p[4]*sVtg + p[5]*p_nT/p_vT;
	
	
	P_T = m_T*cp_AL*T_vT*(1.0 -pow(p_nT/p_vT,(k_AL-1.0)/k_AL))*w_t;
	P_V = m_V*cp_FL*T_vV*(pow(p_nV/p_vV,(k_FL-1.0)/k_FL)-1.0 )*w_v;
	
	
	dx[0] = (P_T - P_V)/(4.0*M_PI*M_PI*J*x_skal)/skal;

	}


	// Optional
	bool ode_structure(DiffStructure &s) override {
	  
//	  return false;

		s(0,x_indexode(0));
		
		s(0,p_indexode(0));
		s(0,p_indexode(1));
		s(0,p_indexode(2));
		s(0,p_indexode(3));
		s(0,p_indexode(4));
		s(0,p_indexode(5));

		return true;

	}

	bool ode_diff(DiffStructure& s, double t, const double* x, const double* u, const double* p) override {

	    double w_t,w_v,P_V,P_T,m_V,T_vV,p_nV, p_vV,m_T,T_vT,p_nT,p_vT,sVtg;  
	    double x_skal = x[0]*skal;
	    
	    //Messdaten auswerten
	    getParams(t,Parameter);
	    
	    m_V = Parameter[0];
	    T_vV = Parameter[1];
	    p_nV = Parameter[2];
	    p_vV = Parameter[3];
	    m_T = Parameter[4];
	    T_vT = Parameter[5];
	    p_nT = Parameter[6];
	    p_vT = Parameter[7];
	    sVtg = Parameter[8];
		    
	    
	    //Auswertung der Polynome
	    w_v = p[0] + p[1]*x[0] + p[2]*p_vV/p_nV;
	    w_t = p[3] + p[4]*sVtg + p[5]*p_nT/p_vT;
	    
	    P_V = m_V*cp_FL*T_vV*(pow(p_nV/p_vV,(k_FL-1.0)/k_FL)-1.0 )*w_v;
	    P_T = m_T*cp_AL*T_vT*(1.0 -pow(p_nT/p_vT,(k_AL-1.0)/k_AL))*w_t;
	    
	    
	    s(0,x_indexode(0)) = (     -(p[1])/(4.0*M_PI*M_PI*J*x_skal) - (P_T - P_V)*(4.0*M_PI*M_PI*J*skal) /((4.0*M_PI*M_PI*J*x_skal)*(4.0*M_PI*M_PI*J*x_skal))    )/skal;

		return true;
	}

	bool ode_diff_p(DiffStructure& s, double t, const double *x, const double *u, const double *p, int index) override {
	    double w_t,w_v,P_V,P_T,m_V,T_vV,p_nV, p_vV,m_T,T_vT,p_nT,p_vT,sVtg;  
	      
	    //return false;
	    double x_skal = x[0]*skal;
	    
	    //Messdaten auswerten
	    getParams(t,Parameter);
	    
	    m_V = Parameter[0];
	    T_vV = Parameter[1];
	    p_nV = Parameter[2];
	    p_vV = Parameter[3];
	    m_T = Parameter[4];
	    T_vT = Parameter[5];
	    p_nT = Parameter[6];
	    p_vT = Parameter[7];
	    sVtg = Parameter[8];
		    
	    
	    //Auswertung der Polynome
	    w_v = p[0] + p[1]*x[0] + p[2]*p_vV/p_nV;
	    w_t = p[3] + p[4]*sVtg + p[5]*p_nT/p_vT;
	    
	    P_T = m_T*cp_AL*T_vT*(1.0 -pow(p_nT/p_vT,(k_AL-1.0)/k_AL))*w_t;
	    P_V = m_V*cp_FL*T_vV*(pow(p_nV/p_vV,(k_FL-1.0)/k_FL)-1.0 )*w_v;
	    
	    
	    s(0,p_indexode(0)) = -(m_V*cp_FL*T_vV*(pow(p_nV/p_vV,(k_FL-1.0)/k_FL)-1.0 )*(1.0))/(4.0*M_PI*M_PI*J*x_skal)/skal;
	    s(0,p_indexode(1)) = -(m_V*cp_FL*T_vV*(pow(p_nV/p_vV,(k_FL-1.0)/k_FL)-1.0 )*(x[0]))/(4.0*M_PI*M_PI*J*x_skal)/skal;;
	    s(0,p_indexode(2)) = -(m_V*cp_FL*T_vV*(pow(p_nV/p_vV,(k_FL-1.0)/k_FL)-1.0 )*(p_vV/p_nV))/(4.0*M_PI*M_PI*J*x_skal)/skal;;
	    s(0,p_indexode(3)) = (m_T*cp_AL*T_vT*(1.0 -pow(p_nT/p_vT,(k_AL-1.0)/k_AL)))/(4.0*M_PI*M_PI*J*x_skal)/skal;
	    s(0,p_indexode(4)) = (m_T*cp_AL*T_vT*(1.0 -pow(p_nT/p_vT,(k_AL-1.0)/k_AL))*(sVtg))/(4.0*M_PI*M_PI*J*x_skal)/skal;
	    s(0,p_indexode(5)) = (m_T*cp_AL*T_vT*(1.0 -pow(p_nT/p_vT,(k_AL-1.0)/k_AL))*(p_nT/p_vT))/(4.0*M_PI*M_PI*J*x_skal)/skal;
		
		return true;
		
	}



	void x_boundary(double *x_low, double *x_upp) override {

		x_low[0] = 1.0E-10;
		x_upp[0] = 1.0E20;



	}

	void p_boundary(double *p_low, double *p_upp) override {
		p_low[0] = -1.0E20;
		p_upp[0] = 1.0E20;
		p_low[1] = -1.0E20;
		p_upp[1] = 1.0E20;
		p_low[2] = -1.0E20;
		p_upp[2] = 1.0E20;
		p_low[3] = -1.0E20;
		p_upp[3] = 1.0E20;
		p_low[4] = -1.0E20;
		p_upp[4] = 1.0E20;
		p_low[5] = -1.0E20;
		p_upp[5] = 1.0E20;
	}
	
	void var_boundary(double *x_low, double *x_upp) override {



		x_low[x_index(0,0)] = 0.4259;
		x_upp[x_index(0,0)] = 0.4259;

	}
	
	void getParams(double t,double *para){
	  

		
		//std::cout<<t_tmp<<endl;

		int idx_u = floor( (t)/dT);
		int idx_o = idx_u+1;
		double t_u = idx_u*dT;
		double Tf = dT*(n_dis-1);
		
	
		 // Alternative lineare Interplation -> kann zu numerischen Probleme führen, da nicht diff'bar
		if (t < Tf){

		  para[0] = (Data.m_vrd[idx_o]-Data.m_vrd[idx_u])/dT*(t-t_u) + Data.m_vrd[idx_u];
		  para[1] = (Data.T_vVrd[idx_o]-Data.T_vVrd[idx_u])/dT*(t-t_u) + Data.T_vVrd[idx_u];
		  para[2] = (Data.p_nVrd[idx_o]-Data.p_nVrd[idx_u])/dT*(t-t_u) + Data.p_nVrd[idx_u];
		  para[3] = (Data.p_vVrd[idx_o]-Data.p_vVrd[idx_u])/dT*(t-t_u) + Data.p_vVrd[idx_u];
		  para[4] = (Data.m_Trb[idx_o]-Data.m_Trb[idx_u])/dT*(t-t_u) + Data.m_Trb[idx_u];
		  para[5] = (Data.T_vTrb[idx_o]-Data.T_vTrb[idx_u])/dT*(t-t_u) + Data.T_vTrb[idx_u];
		  para[6] = (Data.p_nTrb[idx_o]-Data.p_nTrb[idx_u])/dT*(t-t_u) + Data.p_nTrb[idx_u];
		  para[7] = (Data.p_vTrb[idx_o]-Data.p_vTrb[idx_u])/dT*(t-t_u) + Data.p_vTrb[idx_u];
		  para[8] = (Data.sVtg[idx_o]-Data.sVtg[idx_u])/dT*(t-t_u) + Data.sVtg[idx_u];
		  
		}else{
		  para[0] = Data.m_vrd[idx_u];
		  para[1] = Data.T_vVrd[idx_u];
		  para[2] = Data.p_nVrd[idx_u];
		  para[3] = Data.p_vVrd[idx_u];
		  para[4] = Data.m_Trb[idx_u];
		  para[5] = Data.T_vTrb[idx_u];
		  para[6] = Data.p_nTrb[idx_u];
		  para[7] = Data.p_vTrb[idx_u];
		  para[8] = Data.sVtg[idx_u];
		}


	}
	

};


/////////////////////////////////////////////////////////////////////////////



int main(int argv, char* argc[]) {

	TWparameter twparameter("transworhp_Turbolader.xml");
	twparameter.Arguments(argv,argc);
	
		// Force explicit method... 
	#ifdef _USE_EXPL_TransWORHP
	twparameter.twdiscretization = TWdiscretization(TW_MultipleShooting,0,0);
	twparameter.USERHM=-1;
	twparameter.USERDF=-1;
	twparameter.USERDG=-1;
	#endif
	int NDIS = 1001;
	
	double num_node,step;
	if (argv > 1){
	  num_node = ToDouble(argc[1]);
	  if (argv > 2) NDIS = ToDouble(argc[2]);
	  step = NDIS/num_node;
	}else{
	 step = NDIS; 
	}
	
	vector<int> multinode;
	for(int i=0;i<NDIS;i=i+step){
	  multinode.push_back(i);
	}
	
	
	TWfolder folder(&twparameter,0);



	

	Turbo_Expl ph(NDIS,multinode);
	#ifdef _USE_EXPL_TransWORHP
	ph.butcherInit(twparameter.butchertableau,twparameter.stepsize );
	#endif
	ph.LinearTimeAxis(0.0, (NDIS-1)*dT);
	folder.Add(&ph);

	Viewer *viewer = nullptr;
	if (twparameter.PLOT) viewer = new Viewer(&twparameter);

	folder.Init();
	folder.Init(viewer);


	
	//ph.Integrate(twparameter.butchertableau);
	// Ende Optional
	
	folder.Loop(2);
	
	int k;
	std::cout<<"Optimale Parameter:"<<endl;
	for(int i=0;i<6;i++){
	  k = ph.p_index(i);
	  std::cout<<"p1("<<i+1<<") = "<<folder.worhp_o.X[k]<<";"<<endl;
	}
	std::cout<<"Zielfunktionswert = "<<folder.worhp_o.F<<endl;
	
	
	delete viewer;

	return 0;
}

