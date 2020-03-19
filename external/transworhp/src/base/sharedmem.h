#ifndef sharedmem
#define sharedmem

#include "defines.h"
#include <vector>
#include <string>

#ifdef WIN32

#include "windows.h"
extern "C" int writetext_(char *message, int len_message);

#else

int writetext_(char *message, int len_message);

#endif


class SharedMemory
{
	public:
		SharedMemory(int size);
		~SharedMemory();

	
		double &double_sel (int index);
		int &int_sel (int index);
protected:
		int page_size;

#ifdef WIN32
		HANDLE hMapFile;
		char **data;
		int pages;
#else
		int pages;
		int *shmid;
		char **data;
#endif

};

#endif
