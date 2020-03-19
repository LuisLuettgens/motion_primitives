/*----------------------------------------------------------------
 *
 *  Example: Industrieroboter
 *
 *----------------------------------------------------------------*/

#ifdef WIN32
#include "windows.h"
#endif

#ifndef M_PI
#define M_PI 3.14152926535897932846
#endif

#include "TransWORHP.h"

using namespace std;
//using namespace TransWORHP;

vector<double> start;
vector<double> final;


//optctrl -userdg -userhm -userdf +p -n 61 -m 15

TransWorhp *ph=0;

double roboterplot (int &len, int &ndgl, int &nsteuer, double *t,double *x, double *u, int &i, int &index) {

 if (index==0) return x[i*ndgl + 4]-x[i*ndgl + 5];

	if (index==1) {
		int g_index = (ph->n_dis-1) * ph->n_ode + ph->n_rand + i*ph->n_neben;
		return ph->G[g_index + 0];
	}
}

class RoboterPhase : public TransWorhp {
public:
	static string Text[];
	int mode;

	RoboterPhase(int dis, int mode) : TransWorhp(Text[mode].c_str(), dis, 10, 3, 1, 0, 1), mode(mode) {

		freetime = 1;

	}



	void OpenWindows(Viewer *viewer) {

		char buf[30];
		sprintf(buf,"%c_4 - %c_5",130,130);
 viewer->Data(buf, roboterplot,0);
//Funktionenzeiger2 f = 
//		gr->Data(buf, f,-1,1,1);


//		gr->ThreeD("3D");

	}

	void GetXTitle(int d, char *buf) {
		if (d==0) sprintf(buf,"d%c_1",130);
		if (d==1) sprintf(buf,"d%c_2",130);
		if (d==2) sprintf(buf,"d%c_3",130);
		if (d==3) sprintf(buf,"%c_1",130);
		if (d==4) sprintf(buf,"%c_2",130);
		if (d==5) sprintf(buf,"%c_3",130);

		if (d==6) sprintf(buf,"u^2 x^2");
		if (d==7) sprintf(buf,"ddx^2 x^2");
		if (d==8) sprintf(buf,"u^2 + x^2");
		if (d==9) sprintf(buf,"u^2");

	}

	void GetUTitle(int d, char *buf) {
		if (d==0) sprintf(buf,"Control %c_1",130);
		if (d==1) sprintf(buf,"Control %c_2",130);
		if (d==2) sprintf(buf,"Control %c_3",130);

	}

	void p_init(double *p) {

		p[0] = 2.5;

	}

	void init() {

		for (int i=0;i<n_dis;i++) {
			//  X[u_index(i,0)] = 0;
			//  X[u_index(i,1)] = 0;
			//  X[x_index(i,5)] = 5;
			//X[x_index(i,1)] =  + .3 * i/n_dis;
		}

	}

	double obj() {

		//cout << x(n_dis-1,2) << endl;
		if (mode==0) return 0;
		if (mode==1) return x(n_dis-1,6) + .1 * p(0);
		if (mode==2) return x(n_dis-1,7) + .1 * p(0);
		if (mode==3) return x(n_dis-1,8) + .1 * p(0);
		if (mode==4) return x(n_dis-1,9) + .1 * p(0);
		if (mode==5) return p(0);
		return 0;
	}


	bool obj_structure(DiffStructure &s) {

		if (mode==0)
			;

		if (mode==1) {
			s(0,x_index(n_dis-1,6) );
			s(0,p_index(0) );
		}

		if (mode==2) {
			s(0,x_index(n_dis-1,7) );
			s(0,p_index(0) );
		}
		if (mode==3) {
			s(0,x_index(n_dis-1,8) );
			s(0,p_index(0) );
		}
		if (mode==4) {
			s(0,x_index(n_dis-1,9) );
			s(0,p_index(0) );
		}
		if (mode==5) {
			s(0,p_index(0) );
		}

		return true;

	}

	void ode(double *dx, double t_, const double *x, const double *u, const double *p) {

		int ntype = 1;

		double m[] = {510.0,240.0,0};

// Parameter des 1. Arms
		double l1[] = {0.188,0.000,0.900};
		double s1[] = {0.130,0.000,0.783};
		double i1z   =   28.0;

// Parameter des zweiten Arms
		double l2[] = {0.000,0.000,0.950};
		double s2[] = {-0.010,0.007,0.430};
		double i2x  =   33.63;
		double i2y  =   28.43;
		double i2z  =   9.4;

		double WERKZEUG;
		double s3[3];
		double i3x,i3y,i3z;

		if (ntype==1) {

			WERKZEUG = 0.;

			s3[0] =  0.5491;
			s3[1] = -0.0233;
			s3[2] = 0.225 - 0.0079;

			m[2]  =  294.9 + WERKZEUG;
			i3x   =  8.6107;
			i3y   =  181.8004;
			i3z   =  195.2742;
		}
		else if (ntype==2) {

			WERKZEUG = 14.;

			s3[0] =  (294.9*0.5491 + WERKZEUG*(1.7045-0.06))/(294.9+WERKZEUG);
			s3[1] = -0.0233;
			s3[2] =  0.225 - 0.0079;

			m[2]  =  294.9 + WERKZEUG;
			i3x   =  8.6107 + 0.05;
			i3y   =  181.8004;
			i3z   =  195.2742;
		}
		else if (ntype==3) {

			WERKZEUG = 0.;

			s3[0] =  (276.9*0.3697 + WERKZEUG*(1.3-0.06))/(276.9+WERKZEUG);
			s3[1] = -0.0248;
			s3[2] = 0.225 - 0.0106;

			m[2]  =  276.9 + WERKZEUG;
			i3x   =  7.64;
			i3y   =  106.87;
			i3z   =  117.92;
		}

		double sgw[3]={0,0,0};
		double mgw;

		if (ntype==3) {

			sgw[0] = -0.591;
			sgw[2] = -0.041;
			mgw    =  371.0;
		}
		else {

			sgw[0] = -0.665;
			sgw[2] = -0.015;
			mgw    =  465.0;
		}


		double umax[3];
		umax[0] = 16. *36.77 * 14.5 * .4;
		umax[1] = 13. *36.77 * 14.5 * .5;
		umax[2] = 13. *36.77 * 14.5 * .5;

		double UR[3];

		UR[0] = 380. * tanh (3.*x[0]);
		UR[1] = 345. * tanh (3.*x[1]);
		double UH     = 3783.4297*x[4];
		UR[2] = 337. * tanh (3.*x[2]);

		double M[3][3];double rs[3];
		double g = 9.8100;





		double cq3 = cos(x[5]);
		double sq3 = sin(x[5]);

		double cq2 = cos(x[4]);
		double sq2 = sin(x[4]);

		double a1 = s3[0] * cq3 + s3[2] * sq3;
		double a2 = s3[2] * cq3 - s3[0] * sq3;

		double a3 = sgw[0] * cq3 + sgw[2] * sq3;
		double a4 = sgw[2] * cq3 - sgw[0] * sq3;

		double a5 = s2[0] * cq2 + s2[2] * sq2;
		double a6 = s2[2] * cq2 - s2[0] * sq2;

		double  b1 = l2[2] * cq2;
		double b2 = l2[2] * sq2;

		double  c1 = a1*b1-a2*b2;
		double  c2 = a2*b1+a1*b2;

		M[0][0] = m[0]*s1[0]*s1[0]
				  + m[1]*((a5+l1[0])*(a5+l1[0])+s2[1]*s2[1])
				  + m[2]*((a1+b2+l1[0])*(a1+b2+l1[0]) + s3[1]*s3[1]) +
				  mgw*(a3+l1[0])*(a3+l1[0]) + i1z + i2x*sq2*sq2 + i2z*cq2*cq2 + i3x*sq3*sq3 + i3z*cq3*cq3;
		M[0][1] = -m[1]*a6*s2[1] - m[2]*b1*s3[1];
		M[0][2] = -m[2]*a2*s3[1];
		M[1][0] = M[0][1];
		M[1][1] = m[1]*(s2[2]*s2[2] + s2[0]*s2[0]) + m[2]*l2[2]*l2[2] + i2y;
		M[1][2] = m[2]*c2 ;
		M[2][0] = M[0][2];
		M[2][1] = M[1][2];
		M[2][2] = m[2]*(s3[2]*s3[2] + s3[0]*s3[0]) + mgw*(sgw[2]*sgw[2]+sgw[0]*sgw[0])+i3y;

		double d1 = m[1]*(a5+l1[0])*a6 + m[2]*(a1+b2+l1[0])*b1 + (i2x-i2z)*sq2*cq2;
		double d2 = mgw*(a3+l1[0])*a4 + m[2]*(a1+b2+l1[0])*a2 + (i3x-i3z)*sq3*cq3;

		rs[0]  = 2.*x[0]*x[1]*d1 + 2.*x[0] *x[2]*d2 + x[1]*x[1]*(m[1]*a5*s2[1]+m[2]*b2*s3[1])+x[2]*x[2]*m[2]*a1*s3[1];
		rs[1] = -x[0] *x[0]*d1-x[2]*x[2]*m[2]*c1-g*(m[1]*a5+m[2]*b2);
		rs[2] = -x[0] *x[0]*d2+x[1]*x[1]*m[2]*c1-g*(m[2]*a1+mgw*a3);

		rs[0] = -rs[0];
		rs[1] = -rs[1];
		rs[2] = -rs[2];





		rs[0] = rs[0] + u[0]*umax[0] - UR[0];
		rs[1] = rs[1] + u[1]*umax[1] - UR[1] - UH;
		rs[2] = rs[2] + u[2]*umax[2] - UR[2];

		double MM[3][3];
		double t[100];

#include "inv_3dof.c"

		double ddx[3]={0,0,0};

		for (int i=0;i<3;i++) {
			for (int j=0;j<3;j++) {
				ddx[i] += MM[i][j]*rs[j];
			}
		}

		dx[0] = ddx[0];
		dx[1] = ddx[1];
		dx[2] = ddx[2];

		dx[3] = x[0];
		dx[4] = x[1];
		dx[5] = x[2];

		for (int i=6;i<10;i++) {  dx[i] = 0; }

		for (int i=0;i<3;i++) {
			dx[6] +=  u[i] * u[i]  * x[i]*x[i];
			dx[7] +=  ddx[i]*ddx[i] * x[i]*x[i];
			dx[8] +=  x[i]*x[i] + .5*u[i]*u[i];
			dx[9] +=  u[i]*u[i];
		}

		for (int i=0;i<10;i++) {
			dx[i] *= p[0];
		}

	}


// Optional
	bool ode_structure(DiffStructure &s) {

		for (int i=0;i<3;i++) {
			for (int j=0;j<6;j++) {
				s(i,x_index(0,j));
			}
			for (int j=10;j<10+3;j++) {
				s(i,x_index(0,j));
			}
			s(i,p_indexode(0));
		}

		for (int i=3;i<6;i++) {
			s(i,x_index(0,i-3));
			s(i,p_indexode(0));
		}

		// x[6]
		for (int j=0;j<3;j++) {
			s(6,x_index(0,j));
		}
		for (int j=10;j<10+3;j++) {
			s(6,x_index(0,j));
		}
		s(6,p_indexode(0));

		// x[7]
		for (int j=0;j<6;j++) {
			s(7,x_index(0,j));
		}
		for (int j=10;j<10+3;j++) {
			s(7,x_index(0,j));
		}
		s(7,p_indexode(0));

		// x[8]
		for (int j=0;j<3;j++) {
			s(8,x_index(0,j));
		}
		for (int j=10;j<10+3;j++) {
			s(8,x_index(0,j));
		}
		s(8,p_indexode(0));

		// x[9]
		for (int j=10;j<10+3;j++) {
			s(9,x_index(0,j));
		}
		s(9,p_indexode(0));

		return true;
	}

	void u_boundary(double *u_low, double *u_upp) {

		u_low[0] = -1;
		u_upp[0] = +1;

		u_low[1] = -1;
		u_upp[1] = +1;

		u_low[2] = -1;
		u_upp[2] = +1;

		if (mode==0) {
			u_low[0] = 0;
			u_upp[0] = 0;

			u_low[1] = 0;
			u_upp[1] = 0;

			u_low[2] = 0;
			u_upp[2] = 0;
		}
	}

	void x_boundary(double *x_low, double *x_upp) {

		x_low[0] = -100./180.*M_PI;
		x_upp[0] = +100./180.*M_PI;

		x_low[1] = -100./180.*M_PI;
		x_upp[1] = +100./180.*M_PI;

		x_low[2] = -100./180.*M_PI;
		x_upp[2] = +100./180.*M_PI;

		x_low[3] = -1000;
		x_upp[3] = +1000;

		x_low[4] = -70./180.*M_PI;
		x_upp[4] = +70./180.*M_PI;

		x_low[5] = -28./180.*M_PI;
		x_upp[5] = 105./180.*M_PI;

		for (int i=6;i<10;i++) {
			x_low[i] = 0;
	}
	}

	void p_boundary(double *p_low, double *p_upp) {
		p_low[0] = 1;
		p_upp[0] = 20;

		if (mode==0) {
			p_low[0] = 5;
			p_upp[0] = 5;
		}
	}

	void var_boundary(double *x_low, double *x_upp) {

		double start[] = {0,0,0,M_PI/2,0,0};
		double final[] = {0,0,0,-M_PI/2,0,0};


		for (int i=0;i<6;i++) {
			x_low[x_index(0,i)] = x_upp[x_index(0,i)] = start[i];
		}

		for (int i=6;i<10;i++) {
			x_low[x_index(0,i)] = x_upp[x_index(0,i)] = 0;
		}

		if (mode>0) {
			for (int i=0;i<6;i++) {
				x_low[x_index(n_dis-1,i)] = x_upp[x_index(n_dis-1,i)] = final[i];
			}
		}

		x_low[x_index(0,1)] = x_upp[x_index(0,1)] = -1;
		x_low[x_index(n_dis-1,0)] = x_upp[x_index(n_dis-1,0)] = 0;
		x_low[x_index(n_dis-1,1)] = x_upp[x_index(n_dis-1,1)] = 0;
		
	}

	/*void rand(double *r) {
		
		double start[] = {0,0,0,M_PI/2,0,0};
		double final[] = {0,0,0,-M_PI/2,0,0};


		for (int i=0;i<6;i++) {
			r[i] = x(0,i) - start[i];
		}
		for (int i=6;i<10;i++) {
			r[i] = x(0,i);
		}

		for (int i=0;i<6;i++) {
			if (mode==0)
				r[i+10] = 0;
			else
				r[i+10] = x(n_dis-1,i) - final[i];

		}

	}

	*/


	void neben_boundary(double *c_low, double *c_upp) {

		c_low[0] = -.1; //-65./180.*M_PI;
		c_upp[0] = +65./180.*M_PI;

	}

	void neben(double *c, double t, const double *x, const double *u, const double *p) {
		c[0] = x[4]-x[5] + .01*p[0];
			
	}

	bool neben_structure(DiffStructure &s) {

		s(0,x_index(0,4));
		s(0,x_index(0,5));
		s(0,p_indexode(0));

		return true;
	}
	
	bool neben_diff(DiffStructure &s, double t, const double *x, const double *u, const double *p) {

		s(0,x_index(0,4)) =1 ;
		s(0,x_index(0,5)) = -1 ;
		s(0,p_indexode(0))= 0.01;

		return true;
	}
	
	bool neben_diff_p(DiffStructure &s, double t, const double *x, const double *u, const double *p, int index) {
		
		s(0,p_indexode(0)) = .01;
		return true;
	}
	
};

string RoboterPhase::Text[] = {"Roboter Simulation",
							   "Roboter Momente",
							   "Roboter Momente entkoppelt",
							   "Roboter glatt",
							   "Roboter Energie",
							   "Roboter Zeit",
							  };




/////////////////////////////////////////////////////////////////////////////

int main(int argv, char* argc[]) {

	TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv,argc);
	TWfolder folder(&twparameter,0);

	RoboterPhase ph(twparameter.NDIS,1);
	folder.Add(&ph);

	Viewer *viewer = 0;
	if (twparameter.PLOT) viewer = new Viewer(&twparameter);

	folder.Init();
	folder.Init(viewer);
	folder.Loop();

	delete viewer;

	return 0;
}

