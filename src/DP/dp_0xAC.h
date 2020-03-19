#include <iostream>
#include <vector>

#include "udp_client_server.h"
#include "../spline.h"
#include "../global.h"
#include "../utils/functions.h"

#include "dp_base.h"
#include "CourseCorrection.h"

#ifndef dp_0xAC_H
#define dp_0xAC_H


using namespace std;

class dp_0xAC: public dp_base {
	public:
		dp_0xAC();
		void toTW  (double lat0, double lon0, double h0);
};

#endif