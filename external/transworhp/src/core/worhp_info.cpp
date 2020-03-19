#include "worhp_info.h"

#include <iomanip>

using namespace std;

ostream &operator<<(ostream& os, const OptVar &o) {

	for (int i=0;i<o.n;i++) {
		os << "X" << setw(3) << i <<
		setw(10) << o.XL[i] << " <= " <<
		setw(13) << o.X[i] << " <= " <<
		setw(10) << o.XU[i] << endl;
	}
	for (int i=0;i<o.m;i++) {
		os << "G" << setw(3) << i<<
		setw(10) << o.GL[i] << " <= " <<
		setw(13) << o.G[i] << " <= " <<
		setw(10) << o.GU[i] << endl;
	}
	return os;
}

#ifdef WIN32
extern "C" void USERDRAW ( int&,int&,int&,int&,int&,double*,double*,double*,double* ) {}
extern "C" void USERDRAWW ( int&,int&,int&,int&,int&,int&,double*,double*,double*,double*,double* ) {}
extern "C" void USERDRAW2 ( int&,int&,int&,int&,int&,int&,double*,double*,double*,double*,double* ) {}
extern "C" void USERDRAWC ( int&,int&,int&,int&,int&,double*,double*,double*,double* ) {}
#else
extern "C" void userdraw_ ( int&,int&,int&,int&,int&,double*,double*,double*,double* ) {}
extern "C" void userdraww_ ( int&,int&,int&,int&,int&,int&,double*,double*,double*,double*,double* ) {}
extern "C" void userdraw2_ ( int&,int&,int&,int&,int&,int&,double*,double*,double*,double*,double* ) {}
extern "C" void userdrawc_ ( int&,int&,int&,int&,int&,double*,double*,double*,double* ) {}
#endif

void print(WorhpMatrix &M, int index, std::ostream &os) {

	os << M.name << "(";
	if (M.row) os << setw(3) << M.row[index];
	if (M.row && M.col) os << "/" ;
	if (M.col) os << setw(3) << M.col[index];
	os << ") = " ;
	os << setw(10) << M.val[index]<< endl;
}