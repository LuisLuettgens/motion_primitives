#include <iostream>
#include <vector>

#include "udp_client_server.h"
#include "../spline.h"
#include "../global.h"
#include "../utils/functions.h"

#include "dp_base.h"

#ifndef dp_0xA0_H
#define dp_0xA0_H


using namespace std;

class dp_0xA0: public dp_base {
	public:
		dp_0xA0();
		void toTW  (int k, double lat0, double lon0, double h0);
};

#endif