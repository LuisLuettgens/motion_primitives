#include "dp_base.h"

using namespace std;


dp_base::dp_base () {}
dp_base::~dp_base () {}

size_t dp_base::get_size_of_msg() {
	size_t size = 0;
	for(int i=0; i<NBR_ENTRIES_IN_HEADER; i++) size += DTranges[DD[header[i]    .Description].Datatype].Size;
	for(int i=0; i<NBR_ENTRIES_IN_PCKT  ; i++) size += DTranges[DD[packInt[0][i].Description].Datatype].Size * NBR_OF_PCKTS;
	for(int i=0; i<NBR_ENTRIES_IN_FOOTER; i++) size += DTranges[DD[footer[i]    .Description].Datatype].Size;
	
	return size;
}

void dp_base::to_msg( ) {
	char* write_ptr = msg;

	for(int i=0; i<NBR_ENTRIES_IN_HEADER; i++) {
		memcpy(write_ptr, &(header[i].value), DTranges[DD[header[i].Description].Datatype].Size);
		write_ptr += DTranges[DD[header[i].Description].Datatype].Size;
	}
	
	for(int i=0; i<NBR_OF_PCKTS; i++) {
		for(int j=0; j<NBR_ENTRIES_IN_PCKT; j++) {
			memcpy(write_ptr, &(packInt[i][j].value), DTranges[DD[packInt[i][j].Description].Datatype].Size);
			write_ptr += DTranges[DD[packInt[i][j].Description].Datatype].Size;
		}
	}
	
	for(int i=0; i<NBR_ENTRIES_IN_FOOTER; i++) {
		memcpy(write_ptr, &(footer[i].value), DTranges[DD[footer[i].Description].Datatype].Size);
		write_ptr += DTranges[DD[footer[i].Description].Datatype].Size;
	}
}
void dp_base::from_msg ( char* f_msg ) {
	char* read_ptr = f_msg;
	
	for(int i=0; i<NBR_ENTRIES_IN_HEADER; i++) {
		memcpy(&header[i].value, read_ptr,   DTranges[DD[header[i].Description].Datatype].Size);
		read_ptr += DTranges[DD[header[i].Description].Datatype].Size;
		convert(DD[header[i].Description].Datatype, &(header[i].value));
	}
	
	for(int i=0; i<NBR_OF_PCKTS; i++) {
		for(int j=0; j<NBR_ENTRIES_IN_PCKT; j++) {
			memcpy(&packInt[i][j].value, read_ptr,   DTranges[DD[packInt[i][j].Description].Datatype].Size);
			read_ptr += DTranges[DD[packInt[i][j].Description].Datatype].Size;
			convert(DD[packInt[i][j].Description].Datatype, &(packInt[i][j].value));
		}
	}
	
	for(int i=0; i<NBR_ENTRIES_IN_FOOTER; i++) {
		memcpy(&footer[i].value, read_ptr,   DTranges[DD[footer[i].Description].Datatype].Size);
		read_ptr += DTranges[DD[footer[i].Description].Datatype].Size;
		convert(DD[footer[i].Description].Datatype, &(footer[i].value));
	}
	
	toDbl();
}
void dp_base::toInt ( ) {
	for(int i=0; i<NBR_OF_PCKTS; i++) {
		for(int j=0; j<NBR_ENTRIES_IN_PCKT; j++) {
			struct Description d = DD[packInt[i][j].Description];
			double intMin = DTranges[d.Datatype].IntMin;
			double intMax = DTranges[d.Datatype].IntMax;
			packInt[i][j].value = (int)round( intMin + (intMax - intMin)/(d.Max - d.Min)*(pack[i][j] - d.Min) );
		}
	}
}
void dp_base::toDbl ( ) {
	for(int i=0; i<NBR_OF_PCKTS; i++) {
		for(int j=0; j<NBR_ENTRIES_IN_PCKT; j++) {
			struct Description d = DD[packInt[i][j].Description];
			double intMin = DTranges[d.Datatype].IntMin;
			double intMax = DTranges[d.Datatype].IntMax;
			pack[i][j] =  d.Min + (d.Max - d.Min)/(intMax - intMin)*(packInt[i][j].value - intMin);
		}
	}
}

void dp_base::udp_send(const std::string& addr, int port) {
	udp_client_server::udp_client client(addr, port);
	client.send(msg, SIZE_OF_MSG);
}

std::ostream& operator<<(std::ostream& os, const dp_base& datadict){
	ios::fmtflags f( cout.flags() );
	os << endl;
	
	cout << std::setprecision(8);

	for(int i=0; i<datadict.NBR_ENTRIES_IN_HEADER; i++) {
		cout.width(15); os << left << datadict.header[i].Name; cout.width(2); os << left << ": "; cout.width(15); os << right << datadict.header[i].value  << endl;
	}
	for(int i=0; i<datadict.NBR_ENTRIES_IN_FOOTER; i++) {
		cout.width(15); os << left << datadict.footer[i].Name; cout.width(2); os << left << ": "; cout.width(15); os << right << datadict.footer[i].value  << endl;
	}
	os << endl;
	
	cout.width(5); os << right << "Index"; 
	for(int j=0; j<datadict.NBR_ENTRIES_IN_PCKT; j++) {cout.width(15); os << right << datadict.packInt[0][j].Name;}
	os << endl;
	for(int i=0; i<datadict.NBR_OF_PCKTS; i++) {
		cout.width(5); os << right << i; cout << scientific; 
		for(int j=0; j<datadict.NBR_ENTRIES_IN_PCKT; j++) {cout.width(15); os << right << datadict.pack[i][j];}
		os << endl; 
	}
	os << endl;
	
	os << "*********** Encoded Data ***********" << endl;
	cout.width(5); os << right << "Index"; 
	for(int j=0; j<datadict.NBR_ENTRIES_IN_PCKT; j++) {cout.width(15); os << right << datadict.packInt[0][j].Name;}
	os << endl;
	for(int i=0; i<datadict.NBR_OF_PCKTS; i++) {
		cout.width(5); os << right << i; cout << scientific; 
		for(int j=0; j<datadict.NBR_ENTRIES_IN_PCKT; j++) {cout.width(15); os << right << datadict.packInt[i][j].value;}
		os << endl; 
	}
	cout.flags( f );

	return os;
}

void dp_base::fromTW(int k, double lat0, double lon0, double h0, uint32_t initial_time, double delta_time) {

	int n_dis = ph0->n_dis;
	double dtfull; 
	if(ph0->ENDTIME == COMMON_AND_FIXED)       dtfull = ph0->tEnd[0]*ph0->scale_t[0]/(n_dis-1);
	if(ph0->ENDTIME == OPEN_WITH_FIXED_RATIOS) dtfull = ph0->tEnd[0]*ph0->scale_t[0]/(n_dis-1);
	if(ph0->ENDTIME == OPEN)                   dtfull = ph0->tEnd[k]*ph0->scale_t[k]/(n_dis-1);
	double t, dt;
	double low;
	double x0, y0, psi, u, v, r, delta, n;
	
	tk::spline spline_u;
	tk::spline spline_v; 
	tk::spline spline_r;
	
	vector<double> time(n_dis);
	vector<double> uvec(n_dis);
	vector<double> vvec(n_dis);
	vector<double> rvec(n_dis);

	for(int i=0; i<n_dis; i++) {
		if(ph0->ENDTIME == COMMON_AND_FIXED)       t = ph0->tEnd[0]*ph0->scale_t[0]/(n_dis-1)*i; 
		if(ph0->ENDTIME == OPEN_WITH_FIXED_RATIOS) t = ph0->tEnd[0]*ph0->scale_t[0]/(n_dis-1)*i; 
		if(ph0->ENDTIME == OPEN)                   t = ph0->tEnd[k]*ph0->scale_t[k]/(n_dis-1)*i; 
		time[i] = t;
		uvec[i] = ph0->scale[k][3]*ph0->x(i, ph0->NS[k]+3);
		vvec[i] = ph0->scale[k][4]*ph0->x(i, ph0->NS[k]+4);
		rvec[i] = ph0->scale[k][5]*ph0->x(i, ph0->NS[k]+5);
	}
	
	spline_u.set_points(time, uvec);
	spline_v.set_points(time, vvec);
	spline_r.set_points(time, rvec);
	
	t = 0;
	for(int i=0; i<NBR_OF_PCKTS; i++) {
		
		low = (int) floor(t/dtfull);
		dt = t - dtfull*low;
		
		if(low+1 <= n_dis-1) {
			cubicSpline(low, k, dt/dtfull, dtfull, 0, &x0, &y0, &psi, &delta, &n);
		} else {
			cubicSpline(n_dis-1, k, dt, 0, 1, &x0, &y0, &psi, &delta, &n);
		}
		
		u = spline_u.operator()(t);
		v = spline_v.operator()(t);
		r = spline_r.operator()(t);
		
		double lat, lon, h;
		NEDtoLLA(x0, y0, 0, lat0, lon0, h0, &lat, &lon, &h);
		
		pack[i][ 0] = initial_time + t;
		pack[i][ 1] = lat;
		pack[i][ 2] = lon;
		pack[i][ 3] = atan2(sin(psi),cos(psi))>=0 ? atan2(sin(psi),cos(psi))*180./M_PI : (2*M_PI+atan2(sin(psi),cos(psi)))*180./M_PI; 
		pack[i][ 4] = u;
		pack[i][ 5] = v;
		pack[i][ 6] = r*180./M_PI;
		pack[i][ 7] = delta*180./M_PI;
		pack[i][ 8] = 0;
		pack[i][ 9] = n*60;
		pack[i][10] = 0;
		pack[i][11] = 0;
		pack[i][12] = 0;
		pack[i][13] = 0;
		
// 		pack[i][ 0] = initial_time + t;
// 		pack[i][ 1] = lat;
// 		pack[i][ 2] = lon;
// 		pack[i][ 3] = atan2(sin(psi),cos(psi))>=0 ? atan2(sin(psi),cos(psi))*180./M_PI : (2*M_PI+atan2(sin(psi),cos(psi)))*180./M_PI; 
// 		pack[i][ 4] = u;
// 		pack[i][ 5] = v;
// 		pack[i][ 6] = r*180./M_PI;
// 		pack[i][ 7] = delta*180./M_PI;
// 		pack[i][ 8] = n*60;
		
		t += delta_time;
	}
	
	toInt();
	to_msg();
}

void dp_base::toTW  ( ) {}
void dp_base::toTW  (double lat0, double lon0, double h0)  {}
void dp_base::toTW(int k, double lat0, double lon0, double h0)  {}

