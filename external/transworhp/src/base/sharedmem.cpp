#include "sharedmem.h"
#include <iostream>
#include <fstream>
#include <cstring>

using namespace std;

#ifdef WIN32

#include "windows.h"

#ifdef _CONSOLE

int writetext_(char *message, int len_message) {
	MessageBoxA(NULL, message, "TransWORHP", MB_ICONEXCLAMATION);
	return 0;
}

#else

extern "C" int writetext_(char *message, int len_message);

#endif

#else // LINUX

#include <sys/ipc.h>
#include <sys/shm.h>

int writetext_(char *message, int /*len_message*/) {
	cout << "SHAREDMEMORY: " << message << endl;
	return 0;
}

#endif


double d_tmp;
int i_tmp;


SharedMemory::SharedMemory(int size) {

#ifdef WIN32
	pages = 1;
	page_size = size;
#else
	page_size = 1024;
	pages =	(size+1023)/1024;
	shmid = new int[pages];
#endif
	data = new char*[pages];



#ifdef WIN32

	char name[100] = "TRANSWORHP";
	hMapFile = OpenFileMappingA(FILE_MAP_ALL_ACCESS, FALSE, name);

	if (hMapFile == NULL) {
		//writetext_("Could not open file mapping object.",100);// <<            GetLastError() << endl;

		hMapFile = CreateFileMappingA(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, page_size, name);

		if (hMapFile == NULL) {
			writetext_("Could not create or open file mapping object.",100);
			// cout << "Could not create or open file mapping object " << GetLastError() << endl;
			return;
		}
	}


	data[0] =  (char*) MapViewOfFile(hMapFile, FILE_MAP_ALL_ACCESS, 0, 0, page_size);

	if (data[0] == NULL) {
		//   cout << "Could not map view of file " <<            GetLastError() << endl;
		writetext_("Could not map view of file.",100);

		CloseHandle(hMapFile);

		return ;
	}


#else

	key_t *key=new key_t[pages];
	//cout << "Locating shared memory" << endl;
	string id="/tmp/ttt";
	for ( int i=0; i<pages; i++ ) {
		key[i] = ftok ( id.c_str(),'A'+i );
		if ( key[i]==-1 ) {
			throw "Filehandle not found!" ;
			//return; //exit(1);
		}
		shmid[i] = shmget ( key[i],page_size, 0644 | IPC_CREAT );
		if ( shmid[i]==-1 ) cerr << "shmget" <<i <<endl;


		data[i] = ( char* ) shmat ( shmid[i], ( void* ) 0, 0 );
		if ( data[i]== ( char* ) ( -1 ) ) {
			cerr << "shmat" <<i<<endl;
		}

	}

	delete []key;
#endif

}


SharedMemory::~SharedMemory() {

#ifndef WIN32
	for ( int i=0; i<pages; i++ ) {
		int ret = shmdt ( data[i] );
		if ( ret==-1 ) cerr << "shmdt"<<i  <<endl;
	}
	delete []shmid;

#else

	UnmapViewOfFile(data[0]);
	CloseHandle(hMapFile);

#endif

	delete []data;

}


double &SharedMemory::double_sel ( int index ) {

#ifdef WIN32
	if (!hMapFile) return d_tmp;
#endif
	int page = page_size/sizeof ( double );

	int index1 = index%page;
	int index2 = index/page;

	
	if ( index2<pages ) {
		if (data)
			return ( ( double* ) data[index2] ) [index1];
		else return d_tmp;
	}
	
cout << sizeof ( double ) << " " << page << " " <<index  << " " <<index1  << " " <<index2 << endl;

	
	return d_tmp;
}



int &SharedMemory::int_sel ( int index ) {
	
#ifdef WIN32
	if (!hMapFile) return i_tmp;
#endif

	int page = page_size/sizeof ( double );

	int index1 = index%page;
	int index2 = index/page;
	if ( index2<pages ) {
		if (data)
			return ( ( int* ) data[index2] ) [index1];
		else return i_tmp;
	}

	//cout << sizeof ( double ) << " " << page << " " <<index  << " "
	 //    <<index1  << " " <<index2 << endl;

	return i_tmp;
}


