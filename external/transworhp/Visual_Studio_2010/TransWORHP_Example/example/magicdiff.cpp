/*----------------------------------------------------------------
 *
 *
 *
 *----------------------------------------------------------------*/

#ifdef WIN32
#include "windows.h"
#endif

#include "TransWORHP.h"

using namespace std;


/////////////////////////////////////////////////////////////////////////////
#include "newdouble.h"


double desdemonaplot (int &len, int &ndgl, int &nsteuer, double *t,double *x, double *u, int &i, int &index) {


	MagicDouble T(t[i]);
	MagicDouble X(x[i*2]);
	
//	MagicDouble z = T*T;
//	MagicDouble z = sin(T);
//	MagicDouble z = sin(T*T);
//	MagicDouble z = sin(T*T)/(T+.1);
//	MagicDouble z = cos(1/(T+.5));
//	MagicDouble z = exp(-T*T);
	//MagicDouble z = exp(-T*X);
	MagicDouble z = sin(X)+cos(T);
	//cout << z << endl;
	if (index==0) return z.F();
	if (index==1) if (z.DF(&T)) return *z.DF(&T);
	if (index==2) if (z.DF(&X)) return *z.DF(&X);
	
	
	if (index==3) if (z.DDF(&T,&T)) return *z.DDF(&T,&T);
	if (index==4) if (z.DDF(&X,&T)) return *z.DDF(&X,&T);
	if (index==5) if (z.DDF(&T,&X)) return *z.DDF(&T,&X);
	if (index==6) if (z.DDF(&X,&X)) return *z.DDF(&X,&X);
	
	
	return 0;

}

class DummyPhase : public TransWorhp {
public:

	DummyPhase(int dis) : TransWorhp("Dummy",dis,1,1,0,0,0) {
	
		for (int i=0;i<dis;i++) T[i] = (2*M_PI*i)/dis;
		
	}
	
	void OpenWindows(Viewer *viewer) {

		char buf[30];
		sprintf(buf,"f");
		viewer->Data(buf, desdemonaplot,0);
		sprintf(buf,"df dt");
		viewer->Data(buf, desdemonaplot,1);
		sprintf(buf,"df dx");
		viewer->Data(buf, desdemonaplot,2);
		

		sprintf(buf,"d^2f/dt^2");
		viewer->Data(buf, desdemonaplot,3);
		sprintf(buf,"d^2f/dxdt");
		viewer->Data(buf, desdemonaplot,4);
		sprintf(buf,"d^2f/dtdx");
		viewer->Data(buf, desdemonaplot,5);
		sprintf(buf,"d^2f/dx^2");
		viewer->Data(buf, desdemonaplot,6);

	}
	
	double obj() {
		return x(n_dis-1,0);
	}

	void ode(double *dx, double t, const double *x, const double *u, const double *p) {
		dx[0] = 2;
	}

	void u_boundary(double *u_low, double *u_upp) {
		u_low[0] = -1;
		u_upp[0] = +1;
	}
	
	void var_boundary(double *x_low, double *x_upp) {
		x_low[x_index ( 0,0 ) ] = x_upp[x_index ( 0,0 ) ] = -1;
		//x_low[x_index ( n_dis-1,0 ) ] = x_upp[x_index ( n_dis-1,0 ) ] = 1;
	}
};


int main(int argv, char* argc[]) {

	if (1) {
	MagicDouble a(2);
	MagicDouble b(5);
	
	cout << a << endl;
	cout << b << endl;
	cout << "++++" << endl;
	cout << a + b << endl;
	cout << b + a + b + b + b << endl;
	cout << a + 3 << endl;
	cout << 3 + a << endl;
	
	cout << "****" << endl;
	cout << a * b << endl;
	cout << (a + b) * (a+b) << endl;
	cout << a * 3 << endl;
	cout << 3 * a << endl;
	cout << "----" << endl;
	cout << a - b << endl;
//	cout << b + a + b + b + b << endl;
	cout << a - 3 << endl;
	cout << 3 - a << endl;
	cout <<  - a << endl;
	
	cout << "////" << endl;
	cout << a / b << endl;
	cout << (a + b) / (a+b) << endl;
	cout << (a + b) / (a-b) << endl;
	cout << a / 3 << endl;
	cout << 3 / a << endl;
	
	cout<< "::::::::::::::::::" << endl;
	cout << a + b << endl;
	a+=b;
	cout << a << endl;
	

	cout << ":::::::::::::::::::::::::::::::::::" << endl;
	MagicDouble c = a+b*b;
	cout << c << endl;
	MagicDouble d(0);

	cout << c.F() << endl;
	cout << *c.DF(&a) << endl;
	cout << *c.DF(&b) << endl;
	if (c.DF(&d)) cout << *c.DF(&d) << endl;


	}
	else {

	TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv,argc);
	TWfolder folder(&twparameter,0);

	DummyPhase ph(101); 
	folder.Add(&ph);

	Viewer *viewer = 0;
	if (twparameter.PLOT) viewer = new Viewer(&twparameter);

	folder.Init();
	folder.Init(viewer);

	folder.Loop();
	delete viewer;
	
	}

	return 0;
}

