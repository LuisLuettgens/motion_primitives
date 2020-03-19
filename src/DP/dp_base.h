#include <iostream>
#include <vector>

#include "udp_client_server.h"
#include "../spline.h"
#include "../global.h"
#include "../utils/functions.h"

#ifndef dp_base_H
#define dp_base_H


using namespace std;

class dp_base {
	public:
		
		size_t NBR_ENTRIES_IN_HEADER;
		size_t NBR_ENTRIES_IN_PCKT;
		size_t NBR_ENTRIES_IN_FOOTER;
		
		size_t NBR_OF_PCKTS;
		size_t SIZE_OF_MSG;
		
		vector<item> header;
		vector<vector<item>> packInt;
		vector<vector<double>> pack;
		vector<item> footer;
		
		char *msg;

		dp_base();
		~dp_base();
		size_t get_size_of_msg();
		void to_msg   ( );
		void from_msg ( char* f_msg );
		void toInt    ( );
		void toDbl    ( );
		void fromTW(int k, double lat0, double lon0, double h0, uint32_t initial_time, double delta_time);
		void udp_send(const std::string& addr, int port);
		virtual void toTW  ( );
		virtual void toTW  (double lat0, double lon0, double h0);
		virtual void toTW(int k, double lat0, double lon0, double h0);
		friend std::ostream& operator<<(std::ostream& os, const dp_base& datadict);
};

#endif






