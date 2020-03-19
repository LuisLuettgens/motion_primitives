#include <iostream>
#include <vector>

#include "udp_client_server.h"
#include "../spline.h"
#include "../global.h"
#include "../utils/functions.h"

#include "dp_base.h"

#ifndef dp_0xAE_H
#define dp_0xAE_H


using namespace std;

class dp_0xAE: public dp_base {
	public:
		dp_0xAE();
		void toTW  ( );
};

#endif