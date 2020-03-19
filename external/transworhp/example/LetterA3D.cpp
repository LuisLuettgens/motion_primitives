/*----------------------------------------------------------------
 *
 *  Example: Letter A
 *  Author: Sylvain Roy
 *
 *----------------------------------------------------------------*/

#include "TWfolder.h"
#include "TWproblem.h"
#include "TransWORHP.h"

#define n 7
double px[n+1] = {0,0,0,0,0,0,0,0};
double py[n+1] = {0, 11.74, 1.97, 4.61, 11.74, 8.93, 10.61, 15};
double pz[n+1] = {5, 10, 3.22, 0, 10, 2.14, 0, 5};
double theta_init[n] = {60, 180, 270, 360, 80+180, 270, 360};
double theta_final[n] = {0, 270, 360, 360+80, 270, 360, 360+60};
double phi_init[n] = {90,90,90,90,90,90,90};
double phi_final[n] = {90,90,90,90,90,90,90};
double smooth_knot[n-1] = {0,1,1,0,1,1};
int nbr_of_smooth_knots = 4;
double max_v = 3;
int letter_before = 0;
int letter_after = 0;
double a[n] = {1,1,1,1,1,1,1};
double b[n] = {1,1,1,1,1,1,1};
double c[n] = {10,5,1,1,1,1,1};
double d[n] = {1,1,1,1,1,1,1};

using namespace std;

double phasenPlot(int &len, int &node, int &nctrl, double *t,double *x, double *u, int &i, int &index) {
	if (index==0) return x[i*(node+nctrl)+2];
	if (index==1) return x[i*(node+nctrl)+1];
}

class LetterPhase : public tw::TransWorhpProblem {
public:

	int mode;
	
	LetterPhase(int dis, int mode) : TransWorhpProblem(tw::TWdimension("Letter a",dis,12,3,1,0,0)), mode(mode) {
		
	}
	
	string GetXTitle(int d) {
		if (d== 0) return "Position - X";
		if (d== 1) return "Position - Y";
		if (d== 2) return "Position - Z";
		if (d== 3) return "Velocity";
		if (d== 4) return "Angle - Theta";
		if (d== 5) return "Angle - Phi";
		if (d== 6) return "Acceleration";
		if (d== 7) return "Angular Velocity - Theta";
		if (d== 8) return "Angular Velocity - Phi";
		if (d== 9) return "Energy - Velocity";
		if (d==10) return "Energy - Theta";
		if (d==11) return "Energy - Phi";
	}
	
	string GetUTitle(int d) {
		if (d== 0) return "Jerk";
		if (d== 1) return "Angular Acceleration - Theta";
		if (d== 2) return "Angular Acceleration - Phi";
	}
	
	void OpenWindows(tw::Viewer *viewer) override {
		viewer->PhasePlot("Phasen-Plot", phasenPlot,0,1);
	}

	void p_init(double *p) {
		p[0] = 5;
	}	
	
	void x_init(double *x, int i, int dis) {
		for(int k=0; k<n; k++) {
			if(mode == k) {
				  x[0] = px[k] + (px[k+1] - px[k])*i/(dis-1);
				  x[1] = py[k] + (py[k+1] - py[k])*i/(dis-1);
				  x[2] = pz[k] + (pz[k+1] - pz[k])*i/(dis-1);
				  x[3] = 2;
				  x[4] = (theta_init[k] + (theta_final[k] - theta_init[k])*i/(dis-1))*M_PI/180;
				  x[5] = (phi_init[k] + (phi_final[k] - phi_init[k])*i/(dis-1))*M_PI/180;
			}
		}
	}

	double obj() {
		for(int k=0; k<n; k++) {
			if(mode == k) return  a[k]*p(0) + b[k]*x(n_dis-1,9) + c[k]*x(n_dis-1,10) + d[k]*x(n_dis-1,11);
		}
	}

	bool obj_structure(tw::DiffStructure &s) {
		s(0, p_index(0         ));
		s(0, x_index(n_dis-1,9 ));
		s(0, x_index(n_dis-1,10));
		s(0, x_index(n_dis-1,11));
		return true;
	}

	bool obj_diff(tw::DiffStructure &s) {
		for(int k=0; k<n; k++) {
			if(mode == k) {
				s(0, p_index(0         )) = a[k];
				s(0, x_index(n_dis-1,9 )) = b[k];
				s(0, x_index(n_dis-1,10)) = c[k];
				s(0, x_index(n_dis-1,11)) = d[k];
			}
		}
		return true;
	}

	void ode(double *dx, double t, const double *x, const double *u, const double *p) {
		dx[ 0] = x[3]*cos(x[5]);
		dx[ 1] = x[3]*sin(x[5])*cos(x[4]);
		dx[ 2] = x[3]*sin(x[5])*sin(x[4]);
		dx[ 3] = x[6];
		dx[ 4] = x[7];
		dx[ 5] = x[8];
		dx[ 6] = u[0];
		dx[ 7] = u[1];
		dx[ 8] = u[2];
		dx[ 9] = u[0]*u[0];
		dx[10] = u[1]*u[1];
		dx[11] = u[2]*u[2];
	
		for(int i=0; i<12; i++) dx[i] *= p[0];
	}

	bool ode_structure(tw::DiffStructure &s) {
		s( 0, x_indexode(3)); 
		s( 0, x_indexode(5)); 
		s( 1, x_indexode(3)); 
		s( 1, x_indexode(4)); 
		s( 1, x_indexode(5)); 
		s( 2, x_indexode(3)); 
		s( 2, x_indexode(4)); 
		s( 2, x_indexode(5)); 
		s( 3, x_indexode(6));
		s( 4, x_indexode(7));
		s( 5, x_indexode(8));
		s( 6, u_indexode(0));
		s( 7, u_indexode(1));
		s( 8, u_indexode(2));
		s( 9, u_indexode(0));
		s(10, u_indexode(1));
		s(11, u_indexode(2));
		
		for(int i=0; i<12; i++) s(i, p_indexode(0)); 
		return true;
	}

	bool ode_diff(tw::DiffStructure &s, double t, const double *x, const double *u, const double *p) {
		s( 0, x_indexode(3)) = cos(x[5])*p[0]; 
		s( 0, x_indexode(5)) = -x[3]*sin(x[5])*p[0]; 
		s( 1, x_indexode(3)) = sin(x[5])*cos(x[4])*p[0]; 
		s( 1, x_indexode(4)) = -x[3]*sin(x[5])*sin(x[4])*p[0]; 
		s( 1, x_indexode(5)) = x[3]*cos(x[5])*cos(x[4])*p[0]; 
		s( 2, x_indexode(3)) = sin(x[5])*sin(x[4])*p[0]; 
		s( 2, x_indexode(4)) = x[3]*sin(x[5])*cos(x[4])*p[0]; 
		s( 2, x_indexode(5)) = x[3]*cos(x[5])*sin(x[4])*p[0]; 
		s( 3, x_indexode(6)) = 1*p[0];
		s( 4, x_indexode(7)) = 1*p[0];
		s( 5, x_indexode(8)) = 1*p[0];
		s( 6, u_indexode(0)) = 1*p[0];
		s( 7, u_indexode(1)) = 1*p[0];
		s( 8, u_indexode(2)) = 1*p[0];
		s( 9, u_indexode(0)) = 2*u[0]*p[0];
		s(10, u_indexode(1)) = 2*u[1]*p[0];
		s(11, u_indexode(2)) = 2*u[2]*p[0];
		return true;
	}
	
	bool ode_diff_p(tw::DiffStructure &s, double t, const double *x, const double *u, const double *p) {
		s( 0, p_indexode(0)) = x[3]*cos(x[5]);
		s( 1, p_indexode(0)) = x[3]*sin(x[5])*cos(x[4]);
		s( 2, p_indexode(0)) = x[3]*sin(x[5])*sin(x[4]);
		s( 3, p_indexode(0)) = x[6];
		s( 4, p_indexode(0)) = x[7];
		s( 5, p_indexode(0)) = x[8];
		s( 6, p_indexode(0)) = u[0];
		s( 7, p_indexode(0)) = u[1];
		s( 8, p_indexode(0)) = u[2];
		s( 9, p_indexode(0)) = u[0]*u[0];
		s(10, p_indexode(0)) = u[1]*u[1];
		s(11, p_indexode(0)) = u[2]*u[2];
		return true;
	}	
	
	void p_boundary(double *p_low, double *p_upp) {
		p_low[0] = 0.1;
		p_upp[0] = 20;
	}

	void x_boundary(double *x_low, double *x_upp) {
		x_low[ 3] = 0;
		x_upp[ 3] = max_v;
		x_low[ 9] = 0;
		x_low[10] = 0;
		x_low[11] = 0;
	}
	
	void u_boundary(double *u_low, double *u_upp) {
	}

	void var_boundary(double *x_low, double *x_upp) {
	  if(mode == 0) {
	    if(letter_before == 0) {
	      x_low[x_index ( 0,3 ) ] = x_upp[x_index ( 0,3 ) ] = 0;//Initial Velocity       - First Phase
	      x_low[x_index ( 0,6 ) ] = x_upp[x_index ( 0,6 ) ] = 0;//Initial Acceleration   - First Phase
	    } else {
	      x_low[x_index ( 0,3 ) ] = x_upp[x_index ( 0,3 ) ] = max_v;//Initial Velocity   - First Phase
	    }
	    x_low[x_index ( 0,7 ) ] = x_upp[x_index ( 0,7 ) ] = 0;//Initial Angular Velocity - Theta - First Phase
		x_low[x_index ( 0,8 ) ] = x_upp[x_index ( 0,8 ) ] = 0;//Initial Angular Velocity - Phi   - First Phase
	  }
	  if(mode == n-1) {
	    if(letter_after == 0) {
	      x_low[x_index ( n_dis-1,3 ) ] = x_upp[x_index ( n_dis-1,3 ) ] = 0;//Final Velocity       - Last Phase
	      x_low[x_index ( n_dis-1,6 ) ] = x_upp[x_index ( n_dis-1,6 ) ] = 0;//Final Acceleration   - Last Phase
	    } else {
	      x_low[x_index ( n_dis-1,3 ) ] = x_upp[x_index ( n_dis-1,3 ) ] = max_v;//Final Velocity   - Last Phase
	    }
	    x_low[x_index ( n_dis-1,7 ) ] = x_upp[x_index ( n_dis-1,7 ) ] = 0;//Final Angular Velocity - Theta - Last Phase
		x_low[x_index ( n_dis-1,8 ) ] = x_upp[x_index ( n_dis-1,8 ) ] = 0;//Final Angular Velocity - Phi   - Last Phase
	  }

	  for(int k=0; k<n; k++) {
	    if(mode == k) {
	      x_low[x_index ( 0,0 ) ] = x_upp[x_index ( 0,0 ) ] = px[k];//Initial Position - X
	      x_low[x_index ( 0,1 ) ] = x_upp[x_index ( 0,1 ) ] = py[k];//Initial Position - Y
		  x_low[x_index ( 0,2 ) ] = x_upp[x_index ( 0,2 ) ] = pz[k];//Initial Position - Z
	      
	      x_low[x_index ( n_dis-1,0 ) ] = x_upp[x_index ( n_dis-1,0 ) ] = px[k+1];//Final Position - X
	      x_low[x_index ( n_dis-1,1 ) ] = x_upp[x_index ( n_dis-1,1 ) ] = py[k+1];//Final Position - Y
		  x_low[x_index ( n_dis-1,2 ) ] = x_upp[x_index ( n_dis-1,2 ) ] = pz[k+1];//Final Position - Z
	      
	      x_low[x_index ( 0      ,4 ) ] = x_upp[x_index ( 0      ,4 ) ] = theta_init [k]*M_PI/180;//Initial Angle - Theta
		  x_low[x_index ( 0      ,5 ) ] = x_upp[x_index ( 0      ,5 ) ] = phi_init [k]*M_PI/180;//Initial Angle - Phi
	      x_low[x_index ( n_dis-1,4 ) ] = x_upp[x_index ( n_dis-1,4 ) ] = theta_final[k]*M_PI/180;//Final Angle   - Theta
		  x_low[x_index ( n_dis-1,5 ) ] = x_upp[x_index ( n_dis-1,5 ) ] = phi_final[k]*M_PI/180;//Final Angle   - Phi
	      
	      if(k > 0 && smooth_knot[k-1] == 0) {//When previous knot is not smooth
		x_low[x_index ( 0,3 ) ] = x_upp[x_index ( 0,3 ) ] = 0;//Initial Velocity
		x_low[x_index ( 0,6 ) ] = x_upp[x_index ( 0,6 ) ] = 0;//Initial Acceleration
		x_low[x_index ( 0,7 ) ] = x_upp[x_index ( 0,7 ) ] = 0;//Initial Angular Velocity - Theta
		x_low[x_index ( 0,8 ) ] = x_upp[x_index ( 0,8 ) ] = 0;//Initial Angular Velocity - Phi
	      }
	      if(k < n-1 && smooth_knot[k] == 0) {//When next knot is not smooth
		x_low[x_index ( n_dis-1,3 ) ] = x_upp[x_index ( n_dis-1,3 ) ] = 0;//Final Velocity
		x_low[x_index ( n_dis-1,6 ) ] = x_upp[x_index ( n_dis-1,6 ) ] = 0;//Final Acceleration
		x_low[x_index ( n_dis-1,7 ) ] = x_upp[x_index ( n_dis-1,7 ) ] = 0;//Final Angular Velocity - Theta
		x_low[x_index ( n_dis-1,8 ) ] = x_upp[x_index ( n_dis-1,8 ) ] = 0;//Final Angular Velocity - Phi
	      }
	    }
	  }

	  x_low[x_index ( 0, 9 ) ] = x_upp[x_index ( 0, 9 ) ] = 0;//Initial Energy - Velocity 
	  x_low[x_index ( 0,10 ) ] = x_upp[x_index ( 0,10 ) ] = 0;//Initial Energy - Theta
	  x_low[x_index ( 0,11 ) ] = x_upp[x_index ( 0,11 ) ] = 0;//Initial Energy - Phi
	}
	
	void terminate() {
		for(int k=0; k<n; k++) {
			if(mode == k) cout << "End Time of Phase " << k+1 << ": " << p(0) << endl;
		}
	}
};

class LetterAFolder : public tw::TWfolder {
	
public:
	LetterAFolder(tw::TWparameter *p) : TWfolder(p,2*(n-1)+5*(nbr_of_smooth_knots)) {}

	void g_boundary(double *x_low, double *x_upp) {
	
		for(int i=0; i<2*(n-1)+5*(nbr_of_smooth_knots); i++) {
			x_low[i] = 0;
			x_upp[i] = 0;
		}
	}
	
	void con(double *C) {
	  double *X = worhp_o.X;
	  int index1;
	  int index2;
	  int m = 0;
	  
	  for(int k=0; k<n-1; k++) {
	    index1 = phases[k]->x_index(phases[k]->n_dis-1,6) + phases[k  ]->solver->Delta1; //Acceleration
	    index2 = phases[k+1]->x_index(0,6)                + phases[k+1]->solver->Delta1;
	    C[m] = X[index1] - X[index2];
	    
	    index1 = phases[k]->u_index(phases[k]->n_dis-1,0) + phases[k  ]->solver->Delta1; //Jerk
	    index2 = phases[k+1]->u_index(0,0)                + phases[k+1]->solver->Delta1;
	    C[m+1] = X[index1] - X[index2];
	    
	    m = m + 2;
	  }
	  
	  for(int k=0; k<n-1; k++) {
	    if(smooth_knot[k] == 1) {
	      index1 = phases[k]->x_index(phases[k]->n_dis-1,3) + phases[k]->solver->Delta1; //Velocity
	      index2 = phases[k+1]->x_index(0,3)                + phases[k+1]->solver->Delta1;
	      C[m] = X[index1] - X[index2];
	      
	      index1 = phases[k]->x_index(phases[k]->n_dis-1,7) + phases[k]->solver->Delta1; //Angular Velocity - Theta
	      index2 = phases[k+1]->x_index(0,7)                + phases[k+1]->solver->Delta1;
	      C[m+1] = X[index1] - X[index2];
		  
		  index1 = phases[k]->x_index(phases[k]->n_dis-1,8) + phases[k]->solver->Delta1; //Angular Velocity - Phi
	      index2 = phases[k+1]->x_index(0,8)                + phases[k+1]->solver->Delta1;
	      C[m+2] = X[index1] - X[index2];
	      
	      index1 = phases[k]->u_index(phases[k]->n_dis-1,1) + phases[k]->solver->Delta1; //Angular Acceleration - Theta
	      index2 = phases[k+1]->u_index(0,1)                + phases[k+1]->solver->Delta1;
	      C[m+3] = X[index1] - X[index2];
		  
		  index1 = phases[k]->u_index(phases[k]->n_dis-1,2) + phases[k]->solver->Delta1; //Angular Acceleration - Phi
	      index2 = phases[k+1]->u_index(0,2)                + phases[k+1]->solver->Delta1;
	      C[m+4] = X[index1] - X[index2];
	      
	      m = m + 5;
	    }
	  }
	}

	bool con_structure(tw::DiffStructure &s) {
	  int index1;
	  int index2;
	  int m = 0;
	  
	  for(int k=0; k<n-1; k++) {
	    index1 = phases[k]->x_index(phases[k]->n_dis-1,6) + phases[k  ]->solver->Delta1; //Acceleration
	    index2 = phases[k+1]->x_index(0,6)                + phases[k+1]->solver->Delta1;
	    s(m,index1);
	    s(m,index2);
	    
	    index1 = phases[k]->u_index(phases[k]->n_dis-1,0) + phases[k  ]->solver->Delta1; //Jerk
	    index2 = phases[k+1]->u_index(0,0)                + phases[k+1]->solver->Delta1;
	    s(m+1,index1);
	    s(m+1,index2);
	    
	    m = m + 2;
	  }
	  
	  for(int k=0; k<n-1; k++) {
	    if(smooth_knot[k] == 1) {
	      index1 = phases[k]->x_index(phases[k]->n_dis-1,3) + phases[k]->solver->Delta1; // Velocity
	      index2 = phases[k+1]->x_index(0,3)                + phases[k+1]->solver->Delta1;
	      s(m,index1);
	      s(m,index2);
	      
	      index1 = phases[k]->x_index(phases[k]->n_dis-1,7) + phases[k]->solver->Delta1; //Angular Velocity - Theta
	      index2 = phases[k+1]->x_index(0,7)                + phases[k+1]->solver->Delta1;
	      s(m+1,index1);
	      s(m+1,index2);
		  
		  index1 = phases[k]->x_index(phases[k]->n_dis-1,8) + phases[k]->solver->Delta1; //Angular Velocity - Phi
	      index2 = phases[k+1]->x_index(0,8)                + phases[k+1]->solver->Delta1;
	      s(m+2,index1);
	      s(m+2,index2);
	      
	      index1 = phases[k]->u_index(phases[k]->n_dis-1,1) + phases[k]->solver->Delta1; //Angular Acceleration - Theta
	      index2 = phases[k+1]->u_index(0,1)                + phases[k+1]->solver->Delta1;
	      s(m+3,index1);
	      s(m+3,index2);
		  
		  index1 = phases[k]->u_index(phases[k]->n_dis-1,2) + phases[k]->solver->Delta1; //Angular Acceleration - Phi
	      index2 = phases[k+1]->u_index(0,2)                + phases[k+1]->solver->Delta1;
	      s(m+4,index1);
	      s(m+4,index2);
	      
	      m = m + 5;
	    }
	  }
	  return true;
	}

	bool con_diff(tw::DiffStructure &s, int colindex) {
	  int index1;
	  int index2;
	  int m = 0;
	  
	  for(int k=0; k<n-1; k++) {
	    index1 = phases[k]->x_index(phases[k]->n_dis-1,6) + phases[k  ]->solver->Delta1; //Acceleration
	    index2 = phases[k+1]->x_index(0,6)                + phases[k+1]->solver->Delta1;
	    s(m,index1) = +1;
	    s(m,index2) = -1;
	    
	    index1 = phases[k]->u_index(phases[k]->n_dis-1,0) + phases[k  ]->solver->Delta1; //Jerk
	    index2 = phases[k+1]->u_index(0,0)                + phases[k+1]->solver->Delta1;
	    s(m+1,index1) = +1;
	    s(m+1,index2) = -1;
	    
	    m = m + 2;
	  }
	  
	  for(int k=0; k<n-1; k++) {
	    if(smooth_knot[k] == 1) {
	      index1 = phases[k]->x_index(phases[k]->n_dis-1,3) + phases[k]->solver->Delta1; //Velocity
	      index2 = phases[k+1]->x_index(0,3)                + phases[k+1]->solver->Delta1;
	      s(m,index1)= +1;
	      s(m,index2) = -1;
	      
	      index1 = phases[k]->x_index(phases[k]->n_dis-1,7) + phases[k]->solver->Delta1; //Angular Velocity - Theta
	      index2 = phases[k+1]->x_index(0,7)                + phases[k+1]->solver->Delta1;
	      s(m+1,index1)= +1;
	      s(m+1,index2) = -1;
		  
		  index1 = phases[k]->x_index(phases[k]->n_dis-1,8) + phases[k]->solver->Delta1; //Angular Velocity - Phi
	      index2 = phases[k+1]->x_index(0,8)                + phases[k+1]->solver->Delta1;
	      s(m+2,index1)= +1;
	      s(m+2,index2) = -1;
	      
	      index1 = phases[k]->u_index(phases[k]->n_dis-1,1) + phases[k]->solver->Delta1; //Angular Acceleration - Theta
	      index2 = phases[k+1]->u_index(0,1)                + phases[k+1]->solver->Delta1;
	      s(m+3,index1)= +1;
	      s(m+3,index2) = -1;
		  
		  index1 = phases[k]->u_index(phases[k]->n_dis-1,2) + phases[k]->solver->Delta1; //Angular Acceleration - Phi
	      index2 = phases[k+1]->u_index(0,2)                + phases[k+1]->solver->Delta1;
	      s(m+4,index1)= +1;
	      s(m+4,index2) = -1;
	      
	      m = m + 5;
	    }
	  }
	  return true;
	}
	
};

/////////////////////////////////////////////////////////////////////////////

int main(int argv, char* argc[]) {

	tw::TWparameter twparameter("transworhp.xml");
	map<string,string> args = twparameter.Arguments(argv,argc);
	LetterAFolder folder(&twparameter);
	
	LetterPhase *ph[n];
	
	for(int k=0; k<n; k++) {
		ph[k] = new LetterPhase(twparameter.NDIS, k);
		ph[k]->setSolver(&twparameter);
		ph[k]->solver->LinearTimeAxis(k, k+1);
		folder.Add(ph[k]);
	}

	tw::Viewer *viewer = nullptr;
	if (twparameter.PLOT) viewer = new tw::Viewer(&twparameter);

	folder.Init();
	folder.Init(viewer);
	folder.Loop();
	
	char buf[20];
	for(int k=0; k<n; k++) {
		sprintf(buf, "LetterA%d.m", k+1);
		ph[k]->solver->ToMATLAB(buf);
	}
	

	FILE *f = fopen("a11.csv", "w");
	if (f == nullptr) {
		printf("Error opening file!\n");
		exit(1);
	}

	double t;
	double total_t = 0;
	double previous_final_t = 0;
	double vx, vy, vz;
	for(int i=0; i<n-1; i++) {
	  for(int j=0; j<twparameter.NDIS-1; j++) {
	    t = ph[i]->solver->X[ph[i]->p_index(0)]*j/(twparameter.NDIS-1);
	    total_t = previous_final_t + t;
		vx = ph[i]->solver->X[ph[i]->x_index(j,3)]*cos(ph[i]->solver->X[ph[i]->x_index(j,5)]);
	    vy = ph[i]->solver->X[ph[i]->x_index(j,3)]*sin(ph[i]->solver->X[ph[i]->x_index(j,5)])*cos(ph[i]->solver->X[ph[i]->x_index(j,4)]);
	    vz = ph[i]->solver->X[ph[i]->x_index(j,3)]*sin(ph[i]->solver->X[ph[i]->x_index(j,5)])*sin(ph[i]->solver->X[ph[i]->x_index(j,4)]);
	    fprintf(f, "%f, %f, %f, %f\n", total_t, vx, vy, vz);
	  }
	  previous_final_t += ph[i]->solver->X[ph[i]->p_index(0)];
	}
	for(int j=0; j<twparameter.NDIS; j++) {
	    t = ph[n-1]->solver->X[ph[n-1]->p_index(0)]*j/(twparameter.NDIS-1);
	    total_t = previous_final_t + t;
		vx = ph[n-1]->solver->X[ph[n-1]->x_index(j,3)]*cos(ph[n-1]->solver->X[ph[n-1]->x_index(j,5)]);
	    vy = ph[n-1]->solver->X[ph[n-1]->x_index(j,3)]*sin(ph[n-1]->solver->X[ph[n-1]->x_index(j,5)])*cos(ph[n-1]->solver->X[ph[n-1]->x_index(j,4)]);
	    vz = ph[n-1]->solver->X[ph[n-1]->x_index(j,3)]*sin(ph[n-1]->solver->X[ph[n-1]->x_index(j,5)])*sin(ph[n-1]->solver->X[ph[n-1]->x_index(j,4)]);
	    fprintf(f, "%f, %f, %f, %f\n", total_t, vx, vy, vz);
	}
	
	fclose(f);

	for(int k=0; k<n; k++) delete ph[k];
	
	delete viewer;

	return 0;
}
