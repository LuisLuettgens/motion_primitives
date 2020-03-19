#include "vectortools.h"

#ifdef TYPESAFE

int vftablelist[40] = {0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0,
                       0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0
                      };

void addvftable(int a) {

	int i=0;
	while (vftablelist[i]) {
		i++;
	}

	vftablelist[i] = a;


	/*
		i=0;
		while (vftablelist[i]) {
			std::cout << i << " " << (int*)vftablelist[i] << std::endl;
			i++;
		}
	*/

}

#endif
