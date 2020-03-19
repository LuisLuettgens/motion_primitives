#ifndef twsharedmem
#define twsharedmem

#include "../base/sharedmem.h"
#include "TransWORHP.h"

#ifdef WIN32
#include "windows.h"
#endif

/*
const int DATASTART = 1024*12;
const int TEXTPAGE =  8;


class TWsharedmemory : public SharedMemory  {

public:
	TWsharedmemory();

	char *GetText() {
#ifdef WIN32
		if (!hMapFile) return 0;
#endif
		return data[0];	
	}
	char *GetText2(int i) {
#ifdef WIN32
		if (!hMapFile) return 0;
		return &(data[0][0])+i*128;
#else		
		return &data[ i/TEXTPAGE + 1][ (i%TEXTPAGE) * 128];	
#endif
	}
	int &GetIndex() {return int_sel(25);}
	int &GetIteration() {return int_sel(26);}
	int &GetNSTATE() {return int_sel(27);}
	double &GetObj() {return double_sel(28);}
	int &StringIndex() {return int_sel ( 30 );}

	

	void SetData ( TransWorhp *ph);
	void AddString(const std::string &s);
	
	void GetInit ( int &l, int &d1, int &s1, int &n1, int &u1);
	void GetData ( double *t, double *x, double *u, double *g, double *unknown );

private:
	
	int &len() {return int_sel ( 20 );}
	int &ndgl() {return int_sel ( 21 );}
	int &nsteuer() {return int_sel ( 22 );}
	int &nneben() {return int_sel ( 23 );}
	int &nunbe() {return int_sel ( 24 );}

	double &T ( int i )
	{
		return double_sel ( DATASTART +i );
	}
	double &X ( int i, int j )
	{
		return double_sel ( DATASTART + len() + ( ndgl() + nsteuer() ) * i + j );		
	}
	double &U (int i, int j )
	{
		return double_sel ( DATASTART + len() + ( ndgl() + nsteuer() ) * i + ndgl() + j );
		
	}
	double &UNKNOWN ( int i )
	{
		return double_sel ( DATASTART + len() + ( ndgl() + nsteuer() ) * len() + i );
	}

	//int NDISKRET, NSTEUER, NDGL, NUNBE;
	int count;
	int stringpos;
};
*/

#endif

