#include "TWsharedmemory.h"
#include <iostream>
#include <fstream>
#include <cstring>
#include "conversion.h"

using namespace std;

#ifdef WIN32
#include "windows.h"
#else
#include <sys/ipc.h>
#include <sys/shm.h>
#endif

using namespace std;

/*
TWsharedmemory::TWsharedmemory()
	: SharedMemory(1024*100), count(0) {
	
	stringpos = 0;
	
	}


void TWsharedmemory::AddString( const string &message){
	
	vector<string> s = ToStringArray(message, "\n");

	if (s.size()) {
		
		vector<string>::iterator it = s.begin();
	
		
		for (;it!=s.end();it++) {
		//string a = ToString(stringpos) + "::" + it->c_str();
		string a = *it;
		
		memcpy ( GetText2(stringpos), a.c_str(), 128 );
	
		cout << stringpos << "::" << *it << endl;
		
		
		stringpos ++;
		if (stringpos> 8 * 11) stringpos = 0;
		
		StringIndex() = stringpos;
		
		
		
		//MyStatusLine(tag,*it,flag);
		}
	
	
	}


	
}

void TWsharedmemory::SetData(TransWorhp* ph){

	memcpy ( data[0],"TransWORHP",20 );
	//memcpy ( data[1024],"Iteration 1....",20 );
	

	len() = ph->n_dis;
	ndgl() = ph->n_ode;
	nsteuer() = ph->n_ctrl;
	nneben() = ph->n_neben;
	nunbe() = ph->n_param;
	
	for ( int i=0;i< ph->n_dis;i++ ) T ( i ) = ph->T[i];
		
	for ( int i=0;i< ph->n_dis;i++ ) {
		
		for ( int j=0;j< ph->n_ode;j++ ) {
				X ( i,j ) = ph->x(i,j);	
		}
		
		for ( int j=0;j< ph->n_ctrl;j++ ) {
			U ( i,j ) = ph->u(i,j);	
		}
	}
	


	for ( int i=0;i<ph->n_param;i++ ) UNKNOWN ( i ) = ph->p(i);


	
	//	cout << "size:" << len() << " "<< ndgl() << " "<<nsteuer() << "   " <<
	//		len()*(1+ndgl()+nsteuer())*sizeof(double) <<endl;

	GetIndex() = count;
	count++;
}

*/
