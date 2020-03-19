#ifdef WIN32
#include "windows.h"
#endif
#include "optimal_traj_Nships.h"
#include "gridpathsearch.h"
#include "DP/dp_0xA0.h"
#include "DP/dp_0xA1.h"
#include "DP/dp_0xA2.h"
#include "DP/dp_0xA3.h"
#include "DP/dp_0xA4.h"
#include "DP/dp_0xAC.h"
#include "DP/dp_0xAD.h"
#include "DP/dp_0xAE.h"
#include "DP/dp_0xAF.h"

#include <time.h>

TypeOfVehicle vehicle = SHIP;

using namespace std;

int main(int argv, char* argc[]) {
	
	double lat0;
	double lon0;
	double hgt0;
	
	vector<int> shipIDs;
	
	vector <vector<double>> initial_nu;
	vector <vector<double>> final_nu;
	
	if(vehicle == SHIP) {
		lat0 = 54.17057475;
		lon0 = 12.10074142;
		hgt0 = 0;
		
		shipIDs.push_back(0);
		shipIDs.push_back(1);
		
// 		initial_nu.push_back({    54.183615, 12.091236,  +160.0,   0,  0 ,   0 ,   0 ,  0. ,  0 });
// 		initial_nu.push_back({    54.165546, 12.129062,  +270.0,   0,  0 ,   0 ,   0 ,  0. ,  0 });
// 		
// 		final_nu.push_back({    54.160538, 12.118279,   nan("") ,  0 , 0 ,  0 ,  nan("") ,  nan("") ,  nan("") });
// 		final_nu.push_back({    54.148667, 12.093325,   nan("") ,  0,  0 ,  0 ,  nan("") ,  nan("") ,  nan("") });
		
		initial_nu.push_back({    54.151944, 12.123730,  +345.0,   0,  0 ,   0 ,   0 ,  0. ,  0 });
		initial_nu.push_back({    54.187077, 12.088725,  +160.0,   0,  0 ,   0 ,   0 ,  0. ,  0 });
		
		final_nu.push_back({    54.187077, 12.088725,   nan("") ,  0 , 0 ,  0 ,  nan("") ,  nan("") ,  nan("") });
		final_nu.push_back({    54.151944, 12.123730,   nan("") ,  0,  0 ,  0 ,  nan("") ,  nan("") ,  nan("") });
		
	} else {
		lat0 = 53.105988;
		lon0 =  8.852085;
		hgt0 = 0;
		
		shipIDs.push_back(0);
		
		initial_nu.push_back({  53.105008, 8.852455,  +295.0  ,  0,       0,         0,       0,         0,       0  });
		final_nu  .push_back({  53.106884, 8.852185,   nan(""),  0,  nan(""),  nan(""),  nan(""),  nan(""),  nan("") });
	}
	
	{
		dp_0xAC datadict_0xAC;
		datadict_0xAC.header[0].value = 1;
		datadict_0xAC.header[1].value = DP_0XAC;
		datadict_0xAC.header[2].value = initial_nu.size();
		datadict_0xAC.header[3].value = 1000;
		
		for(int k=0; k<datadict_0xAC.header[2].value; k++) {

			datadict_0xAC.pack[k][ 0] = k;
			datadict_0xAC.pack[k][ 1] = initial_nu[k][0];
			datadict_0xAC.pack[k][ 2] = initial_nu[k][1];
			datadict_0xAC.pack[k][ 3] = initial_nu[k][2];
			datadict_0xAC.pack[k][ 4] = initial_nu[k][3];
			datadict_0xAC.pack[k][ 5] = initial_nu[k][4];
			datadict_0xAC.pack[k][ 6] = initial_nu[k][5];
			
			datadict_0xAC.pack[k][ 7] = initial_nu[k][6];
			datadict_0xAC.pack[k][ 8] = 0;
			datadict_0xAC.pack[k][ 9] = initial_nu[k][7];
			datadict_0xAC.pack[k][10] = 0;
			datadict_0xAC.pack[k][11] = 0;
			datadict_0xAC.pack[k][12] = 0;
			datadict_0xAC.pack[k][13] = 0;

			datadict_0xAC.pack[k][14] = 0;
			datadict_0xAC.pack[k][15] = 0;
			datadict_0xAC.pack[k][16] = 0;
			datadict_0xAC.pack[k][17] = 0;
			
			datadict_0xAC.pack[k][18] = 0;
			
			datadict_0xAC.pack[k][19] = 0;
			datadict_0xAC.pack[k][20] = 0;
			datadict_0xAC.pack[k][21] = 0;
			datadict_0xAC.pack[k][22] = 0;
			datadict_0xAC.pack[k][23] = 0;
			datadict_0xAC.pack[k][24] = 0;
			
			datadict_0xAC.pack[k][25] = 0;
			datadict_0xAC.pack[k][26] = 0;
			datadict_0xAC.pack[k][27] = 0;
			datadict_0xAC.pack[k][28] = 0;
		}
		
		datadict_0xAC.toInt();
		datadict_0xAC.to_msg();
		std::this_thread::sleep_for (std::chrono::milliseconds(100));
		cout << "################## Sending message" << "  " << time(NULL) << endl;
		datadict_0xAC.udp_send(ipad, port2);
		cout << "################## Message sent" << "  " << time(NULL) << endl;
	}
	
	{
		dp_0xAD datadict_0xAD;
		datadict_0xAD.header[0].value = 1;
		datadict_0xAD.header[1].value = DP_0XAD;
		datadict_0xAD.header[2].value = 1;
		datadict_0xAD.header[3].value = 1000;
		
		for(int k=0; k<datadict_0xAD.header[2].value; k++) {
				datadict_0xAD.pack[k][ 0] = k;
				datadict_0xAD.pack[k][ 1] = 200.0;
				datadict_0xAD.pack[k][ 2] =  47.5;
				if(vehicle == SHIP) {
					datadict_0xAD.pack[k][ 3] = 54.148597;
					datadict_0xAD.pack[k][ 4] = 12.111736;
				} else {
					datadict_0xAD.pack[k][ 3] = 53.105988;
					datadict_0xAD.pack[k][ 4] = 8.852085;
				}
				datadict_0xAD.pack[k][ 5] = 345.;
				datadict_0xAD.pack[k][ 6] = 345.;
				datadict_0xAD.pack[k][ 7] = 4.*0.514444;
				datadict_0xAD.pack[k][ 8] = 0;
				datadict_0xAD.pack[k][ 9] = 0;
				datadict_0xAD.pack[k][10] = 0;
		}
		
		datadict_0xAD.toInt();
		datadict_0xAD.to_msg();
		std::this_thread::sleep_for (std::chrono::milliseconds(100));
		cout << "################## Sending message" << "  " << time(NULL) << endl;
		datadict_0xAD.udp_send(ipad, port2);
		cout << "################## Message sent" << "  " << time(NULL) << endl;
	}

	{
		dp_0xAE datadict_0xAE;
		datadict_0xAE.header[0].value = 1;
		datadict_0xAE.header[1].value = DP_0XAE;
		datadict_0xAE.header[2].value = 1;
		datadict_0xAE.header[3].value = 1000;
	
		for(int k=0; k<datadict_0xAE.header[2].value; k++) {

			double lat_;
			double lon_;
			if(vehicle == SHIP) {
				lat_ = 54.169407;
				lon_ = 12.099955;
			} else {
				lat_ = 53.105988;
				lon_ = 8.852085;
			}
			
			datadict_0xAE.pack[k][ 0] = lat_;
			datadict_0xAE.pack[k][ 1] = lon_;
			datadict_0xAE.pack[k][ 2] = 0;
			datadict_0xAE.pack[k][ 3] = 0;
			datadict_0xAE.pack[k][ 4] = 0;
			datadict_0xAE.pack[k][ 5] = 0;
		}
		
		datadict_0xAE.toInt();
		datadict_0xAE.to_msg();
		std::this_thread::sleep_for (std::chrono::milliseconds(100));
		cout << "################## Sending message" << "  " << time(NULL) << endl;
		datadict_0xAE.udp_send(ipad, port2);
		cout << "################## Message sent" << "  " << time(NULL) << endl;
	}

	
	for(int k=0; k<initial_nu.size(); k++) {
		dp_0xAF datadict_0xAF;
		datadict_0xAF.header[0].value = 1;
		datadict_0xAF.header[1].value = DP_0XAF;
		datadict_0xAF.header[2].value = k;
		
		datadict_0xAF.pack[0][ 0] = final_nu[k][0];
		datadict_0xAF.pack[0][ 1] = final_nu[k][1];
		datadict_0xAF.pack[0][ 2] = final_nu[k][2];
		datadict_0xAF.pack[0][ 3] = final_nu[k][3];
		datadict_0xAF.pack[0][ 4] = final_nu[k][4];
		datadict_0xAF.pack[0][ 5] = final_nu[k][5];
		datadict_0xAF.pack[0][ 6] = 3;
		
		datadict_0xAF.toInt();
		datadict_0xAF.to_msg();
		std::this_thread::sleep_for (std::chrono::milliseconds(100));
		cout << "################## Sending message" << "  " << time(NULL) << endl;
		datadict_0xAF.udp_send(ipad, port2);
		cout << "################## Message sent" << "  " << time(NULL) << endl;
	}
	

	double dt = 0;
	while(1) {
		dt += 5;
		
		dp_0xAC datadict_0xAC;
		datadict_0xAC.header[0].value = 1;
		datadict_0xAC.header[1].value = DP_0XAC;
		datadict_0xAC.header[2].value = initial_nu.size();
		datadict_0xAC.header[3].value = 1000;
		for(int k=0; k<datadict_0xAC.header[2].value; k++) {
			uint8_t packetID;
			udp_client_server::udp_server server(ipad, port1);

			char recv_msg[MaxsizeofMsg] = {'\0'};
			cout << "################## Waiting for message" << "  " << time(NULL) << endl;
			int sizeofMsg = server.recv(recv_msg, MaxsizeofMsg);
			cout << "################## Message received" << "  " << time(NULL) << endl;

			memcpy(&packetID, recv_msg+2, 1); 


			if(packetID == DP_0XA3) {
				dp_0xA3 datadict_0xA3;
				datadict_0xA3.from_msg(recv_msg);

				double lat, lon, psi, u, v, r, delta, n;

				tk::spline spline_lat;
				tk::spline spline_lon; 
				tk::spline spline_psi;
				tk::spline spline_u;
				tk::spline spline_v; 
				tk::spline spline_r;
				tk::spline spline_delta; 
				tk::spline spline_n;

				int n_dis = 40;
				vector<double> time    (n_dis);
				vector<double> latvec  (n_dis);
				vector<double> lonvec  (n_dis);
				vector<double> psivec  (n_dis);
				vector<double> uvec    (n_dis);
				vector<double> vvec    (n_dis);
				vector<double> rvec    (n_dis);
				vector<double> deltavec(n_dis);
				vector<double> nvec    (n_dis);

				for(int i=0; i<n_dis; i++) {
					time[i]     = datadict_0xA3.pack[i][0] - datadict_0xA3.pack[0][0];
					latvec[i]   = datadict_0xA3.pack[i][1];
					lonvec[i]   = datadict_0xA3.pack[i][2];
					psivec[i]   = datadict_0xA3.pack[i][3];
					uvec[i]     = datadict_0xA3.pack[i][4];
					vvec[i]     = datadict_0xA3.pack[i][5];
					rvec[i]     = datadict_0xA3.pack[i][6];
					deltavec[i] = datadict_0xA3.pack[i][7];
					nvec[i]     = datadict_0xA3.pack[i][9];
				}

				spline_lat  .set_points(time, latvec);
				spline_lon  .set_points(time, lonvec);
				spline_psi  .set_points(time, psivec);
				spline_u    .set_points(time, uvec);
				spline_v    .set_points(time, vvec);
				spline_r    .set_points(time, rvec);
				spline_delta.set_points(time, deltavec);
				spline_n    .set_points(time, nvec);
				
// 				double t_interv = 5; 
				double t_interv = 0; 
				
				lat    = spline_lat .operator()(t_interv);
				lon    = spline_lon .operator()(t_interv);
				psi   = spline_psi  .operator()(t_interv);
				u     = spline_u    .operator()(t_interv);
				v     = spline_v    .operator()(t_interv);
				r     = spline_r    .operator()(t_interv);
				delta = spline_delta.operator()(t_interv);
				n     = spline_n    .operator()(t_interv);


				datadict_0xAC.pack[k][ 0] = k;
				
				datadict_0xAC.pack[k][ 1] = lat;
				datadict_0xAC.pack[k][ 2] = lon;
				datadict_0xAC.pack[k][ 3] = psi; 
				datadict_0xAC.pack[k][ 4] = u;
				datadict_0xAC.pack[k][ 5] = v;
				datadict_0xAC.pack[k][ 6] = r;
				
				datadict_0xAC.pack[k][ 7] = delta;
				datadict_0xAC.pack[k][ 8] = 0;
				datadict_0xAC.pack[k][ 9] = n;
				datadict_0xAC.pack[k][10] = 0;
				datadict_0xAC.pack[k][11] = 0;
				datadict_0xAC.pack[k][12] = 0;
				datadict_0xAC.pack[k][13] = 0;

				datadict_0xAC.pack[k][14] = 0;
				datadict_0xAC.pack[k][15] = 0;
				datadict_0xAC.pack[k][16] = 0;
				datadict_0xAC.pack[k][17] = 0;
				
				datadict_0xAC.pack[k][18] = 0;
				
				datadict_0xAC.pack[k][19] = 0;
				datadict_0xAC.pack[k][20] = 0;
				datadict_0xAC.pack[k][21] = 0;
				datadict_0xAC.pack[k][22] = 0;
				datadict_0xAC.pack[k][23] = 0;
				datadict_0xAC.pack[k][24] = 0;
				
				datadict_0xAC.pack[k][25] = 0;
				datadict_0xAC.pack[k][26] = 0;
				datadict_0xAC.pack[k][27] = 0;
				datadict_0xAC.pack[k][28] = 0;
			}
			if(packetID == DP_0XA2) {
				dp_0xA2 datadict_0xA2;
				datadict_0xA2.from_msg(recv_msg);

				double lat, lon, psi, u, v, r, delta, n;

				tk::spline spline_lat;
				tk::spline spline_lon; 
				tk::spline spline_psi;
				tk::spline spline_u;
				tk::spline spline_v; 
				tk::spline spline_r;
				tk::spline spline_delta; 
				tk::spline spline_n;

				int n_dis = 40;
				vector<double> time    (n_dis);
				vector<double> latvec  (n_dis);
				vector<double> lonvec  (n_dis);
				vector<double> psivec  (n_dis);
				vector<double> uvec    (n_dis);
				vector<double> vvec    (n_dis);
				vector<double> rvec    (n_dis);
				vector<double> deltavec(n_dis);
				vector<double> nvec    (n_dis);

				for(int i=0; i<n_dis; i++) {
					time[i]     = datadict_0xA2.pack[i][0] - datadict_0xA2.pack[0][0];
					latvec[i]   = datadict_0xA2.pack[i][1];
					lonvec[i]   = datadict_0xA2.pack[i][2];
					psivec[i]   = datadict_0xA2.pack[i][3];
					uvec[i]     = datadict_0xA2.pack[i][4];
					vvec[i]     = datadict_0xA2.pack[i][5];
					rvec[i]     = datadict_0xA2.pack[i][6];
					deltavec[i] = datadict_0xA2.pack[i][7];
					nvec[i]     = datadict_0xA2.pack[i][9];
				}

				spline_lat  .set_points(time, latvec);
				spline_lon  .set_points(time, lonvec);
				spline_psi  .set_points(time, psivec);
				spline_u    .set_points(time, uvec);
				spline_v    .set_points(time, vvec);
				spline_r    .set_points(time, rvec);
				spline_delta.set_points(time, deltavec);
				spline_n    .set_points(time, nvec);
				
				int t_interv = 2; 
				
				
				lat    = spline_lat .operator()(t_interv);
				lon    = spline_lon .operator()(t_interv);
				psi   = spline_psi  .operator()(t_interv);
				u     = spline_u    .operator()(t_interv);
				v     = spline_v    .operator()(t_interv);
				r     = spline_r    .operator()(t_interv);
				delta = spline_delta.operator()(t_interv);
				n     = spline_n    .operator()(t_interv);

				datadict_0xAC.pack[k][ 0] = k;
				
				datadict_0xAC.pack[k][ 1] = lat;
				datadict_0xAC.pack[k][ 2] = lon;
				datadict_0xAC.pack[k][ 3] = psi; 
				datadict_0xAC.pack[k][ 4] = u;
				datadict_0xAC.pack[k][ 5] = v;
				datadict_0xAC.pack[k][ 6] = r;
				
				datadict_0xAC.pack[k][ 7] = delta;
				datadict_0xAC.pack[k][ 8] = 0;
				datadict_0xAC.pack[k][ 9] = n;
				datadict_0xAC.pack[k][10] = 0;
				datadict_0xAC.pack[k][11] = 0;
				datadict_0xAC.pack[k][12] = 0;
				datadict_0xAC.pack[k][13] = 0;

				datadict_0xAC.pack[k][14] = 0;
				datadict_0xAC.pack[k][15] = 0;
				datadict_0xAC.pack[k][16] = 0;
				datadict_0xAC.pack[k][17] = 0;
				
				datadict_0xAC.pack[k][18] = 0;
				
				datadict_0xAC.pack[k][19] = 0;
				datadict_0xAC.pack[k][20] = 0;
				datadict_0xAC.pack[k][21] = 0;
				datadict_0xAC.pack[k][22] = 0;
				datadict_0xAC.pack[k][23] = 0;
				datadict_0xAC.pack[k][24] = 0;
				
				datadict_0xAC.pack[k][25] = 0;
				datadict_0xAC.pack[k][26] = 0;
				datadict_0xAC.pack[k][27] = 0;
				datadict_0xAC.pack[k][28] = 0;
			}
		}
		datadict_0xAC.toInt();
		datadict_0xAC.to_msg();
		std::this_thread::sleep_for (std::chrono::milliseconds(100));
		cout << "################## Sending message" << "  " << time(NULL) << endl;
		datadict_0xAC.udp_send(ipad, port2);
		cout << "################## Message sent" << "  " << time(NULL) << endl;
		
		
		{
			dp_0xAD datadict_0xAD;
			datadict_0xAD.header[0].value = 1;
			datadict_0xAD.header[1].value = DP_0XAD;
			datadict_0xAD.header[2].value = 1;
			datadict_0xAD.header[3].value = 1000;
			
			double x0, y0, z0;
			for(int k=0; k<datadict_0xAD.header[2].value; k++) {
					
					double lat;
					double lon;
					double hgt;
// 					double lat_ = 54.161474;
// 					double lon_ = 12.097900;
					double lat_;
					double lon_;
					if(vehicle == SHIP) {
						lat_ = /*54.161474;*/54.160463;
						lon_ = 12.097900;
					} else {
						lat_ = 53.105988;
						lon_ = 8.852085;
					}
					double psi_ = 90;
					double sog_ = 2.0*0.514444;
					LLAtoNED(lat_, lon_, 0, lat0, lon0, hgt0, &x0, &y0, &z0);
					x0 = x0 + dt*sog_*cos(psi_*M_PI/180);
					y0 = y0 + dt*sog_*sin(psi_*M_PI/180);
					NEDtoLLA(x0, y0, 0, lat0, lon0, hgt0, &lat, &lon, &hgt);
					
					datadict_0xAD.pack[k][ 0] = k;
					datadict_0xAD.pack[k][ 1] = 200.0;
					datadict_0xAD.pack[k][ 2] =  47.5;
					datadict_0xAD.pack[k][ 3] = lat;
					datadict_0xAD.pack[k][ 4] = lon;
					datadict_0xAD.pack[k][ 5] = psi_;
					datadict_0xAD.pack[k][ 6] = psi_;
					datadict_0xAD.pack[k][ 7] = sog_;
					datadict_0xAD.pack[k][ 8] = 0;
					datadict_0xAD.pack[k][ 9] = 0;
					datadict_0xAD.pack[k][10] = 0;
			}
			
			datadict_0xAD.toInt();
			datadict_0xAD.to_msg();
			std::this_thread::sleep_for (std::chrono::milliseconds(100));
			cout << "################## Sending message" << "  " << time(NULL) << endl;
			datadict_0xAD.udp_send(ipad, port2);
			cout << "################## Message sent" << "  " << time(NULL) << endl;
		}
		
		{
			dp_0xAE datadict_0xAE;
			datadict_0xAE.header[0].value = 1;
			datadict_0xAE.header[1].value = DP_0XAE;
			datadict_0xAE.header[2].value = 1;
			datadict_0xAE.header[3].value = 1000;
			
			for(int k=0; k<datadict_0xAE.header[2].value; k++) {

// 					double lat_ = 54.169407;
// 					double lon_ = 12.099955;
					double lat_;
					double lon_;
					if(vehicle == SHIP) {
// 						lat_ = 54.169407;
// 						lon_ = 12.099955;
						lat_ = 54.150586; 
						lon_ = 12.117273;
					} else {
						lat_ = 53.105988;
						lon_ = 8.852085;
					}

					datadict_0xAE.pack[k][ 0] = lat_;
					datadict_0xAE.pack[k][ 1] = lon_;
					datadict_0xAE.pack[k][ 2] = 0;
					datadict_0xAE.pack[k][ 3] = 0;
					datadict_0xAE.pack[k][ 4] = 0;
					datadict_0xAE.pack[k][ 5] = 0;
			}
			
			datadict_0xAE.toInt();
			datadict_0xAE.to_msg();
			std::this_thread::sleep_for (std::chrono::milliseconds(100));
			cout << "################## Sending message" << "  " << time(NULL) << endl;
			datadict_0xAE.udp_send(ipad, port2);
			cout << "################## Message sent" << "  " << time(NULL) << endl;
		}
		
	}
	
	return 0;
}