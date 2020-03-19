#include "xopt_data.h"

#include <stdio.h>
#include <iostream>

namespace tw {

int DataStorage::timedgl = -1;
DataStorage::TimeMode_e DataStorage::timemode = dgl_as_time;
double DataStorage::floattime;
double DataStorage::floatstart;

int DataStorage::nsplitdis = 0;
int DataStorage::nsplitlen = 0;
int DataStorage::nsplitdis2 = 0;
int *DataStorage::nsplitindex = nullptr;
double DataStorage::MaxTF = 0;
double DataStorage::MinT0 = 0;

/*extern "C" void settimeistime_() {

	DataStorage::timemode = DataStorage::time_is_time;

	//std::cout << "settimeistime " << endl;
}*/

extern "C" void settimedgl_(int &timedgl_) {

	DataStorage::timedgl = timedgl_;
	DataStorage::timemode = DataStorage::dgl_as_time;

	//std::cout << "settimedgl " << endl;
}
extern "C" void settimedglsplit_(int &timedgl_, int &splitd, int &splitl, int &splitd2, int &spliti) {

	DataStorage::timedgl = timedgl_;
	DataStorage::timemode = DataStorage::dgl2_as_time;
	DataStorage::nsplitdis = splitd;
	DataStorage::nsplitlen = splitl;
	DataStorage::nsplitdis2 = splitd2;
	DataStorage::nsplitindex = &spliti;
	/*for (int i=0;i<splitl;i++) {
		std::cout <<"SSDSSSSSSSSSS"<< DataStorage::nsplitindex[i] << std::endl;
	}*/
}
extern "C" void setfixedtime_() {

	DataStorage::timemode = DataStorage::index_as_time;

	//std::cout << "setfixedtime " << endl;
}


extern "C" void settimefunc_() {

	DataStorage::timemode = DataStorage::time_as_func;

	//TimeFunc = 0;
	//timefunc_(i);
	//std::cout << "settimefunc " << endl;
}

timefuncf TimeFunc = nullptr;

extern "C" void settimefunc2_(timefuncf f) {

	TimeFunc = f;
	//DataStorage::timemode = DataStorage::time_as_func;

	//TimeFunc = 0;
	//timefunc_(i);
	//std::cout << "settimefunc " << endl;
}



double CallTimeFunc(int &i,double *X, double *UNBE, double *T, int &len, int &ndgl, int &nunbe) {
	
	if (TimeFunc) {
		return TimeFunc(i,X,UNBE,T,len,ndgl,nunbe);
	} else {
		return 0;
	}
}

/*extern "C" void setfloattime_(double &t) {

	DataStorage::timemode = DataStorage::fixed_float_as_time;
	DataStorage::floattime=t;
//	std::cout << "setfloattime " << t << endl;
}*/


DataVector::DataVector() : ndis(0), size(0), offset(0), D(nullptr) { }
	

	
void DataVector::setVector(int n, int s, int off, double *value) {
	ndis = n;
	size = s;
	offset = off;
	D = value;
}
double DataVector::getData(int dgl, int index) const {
	return D[offset*index+dgl];
}

void DataVector::setData(int dgl, int index, double v) {
	if (dgl < size && dgl >= 0 && index < ndis && index >= 0) {
		D[offset*index+dgl] = v;
		//std::cout << "D: " << dgl << " " << index << " " << v << std::endl;
	}
}
	
void DataVector::freeze() {
	auto tmp = new double[ndis*offset];
	for (int i = 0; i < ndis*offset; i++) {
		tmp[i] = D[i];
	}
	D = tmp;
}

void DataVector::unFreeze() {
	delete D;
	D = nullptr;
	ndis = 0;
	size = 0;
}

void DataVector::copy(DataVector &other) {
	ndis = other.ndis;
	size = other.size;
	offset = other.offset;
	for (int i = 0; i < ndis*offset; i++) {
		D[i] = other.D[i];
	}
}

int DataVector::getSize() const {
	return size;
}



DataStorage::DataStorage() : tfgesamt(1.0), t0gesamt(0.0), length(0), disLengthX(0), disLengthU(0) {}

	
void DataStorage::setData(int len, double *t,
			  double *x, int xoffset, int xdata,
			  double *u, int uoffset, int udata,
			  double *g, int goffset, int gdata,
			  double *l, int loffset, int ldata) {
	
	length = len;

	T.setVector(len,1,1,t);

	X.setVector(len,xdata,xoffset,x);
	U.setVector(len,udata,uoffset,u);
	G.setVector(len,gdata,goffset,g);
	L.setVector(len,ldata,loffset,l);

	//tfgesamt = X.GetData(size-1,ndgl-1);

	switch (timemode) {
		case fixed_float_as_time:
			tfgesamt = floattime;
			start = floatstart;
			break;

		case index_as_time:
			tfgesamt = len;
			break;

		case time_is_time:
			tfgesamt = T.getData(0,len-1);
			t0gesamt = T.getData(0,0);
			break;

		case time_as_func: {
			int nn = len - 1;
			tfgesamt = CallTimeFunc(nn,X.D,UNKNOWN.D,T.D,X.ndis,X.size,UNKNOWN.size);
			}
			break;

		case dgl_as_time:
			if (timedgl==-1) {
				tfgesamt = x[len*xoffset-1];
			} else {
				tfgesamt = X.getData(timedgl-1,1);
			}
			break;

		case dgl2_as_time:
			//if (timedgl==-1) {
			//	tfgesamt = X.GetData(X.GetSize()-1,1);
			//} else
		{
			tfgesamt = X.getData(timedgl-1,1);
		}
		break;
	}

	//	cout << "TFGESAMT = " << tfgesamt << endl;
}


void DataStorage::addDisData(int disLenX, int disLenU, double *dis_t, double *dis_x, double *dis_u) {

	disLengthX = disLenX;
	disLengthU = disLenU;

	disT.setVector(disLenX,1,1,dis_t);
	disX.setVector(disLenX,X.size,X.offset,dis_x);
	disU.setVector(disLenU,U.size,U.offset,dis_u);
}


void DataStorage::addDisData(int disLenX, int disLenU, int xoffset, int uoffset, double *dis_t, double *dis_x, double *dis_u) {

	disLengthX = disLenX;
	disLengthU = disLenU;

	disT.setVector(disLenX,1,1,dis_t);
	disX.setVector(disLenX,X.size,xoffset,dis_x);
	disU.setVector(disLenU,U.size,uoffset,dis_u);
}


void DataStorage::setData(const Selector &data, int index, double val) {

	switch (data.c) {
		case 't':
			T.setData(data.n,index,val);
			break;
		case 'x':
			X.setData(data.n,index,val);
			break;
		case 'u':
			U.setData(data.n,index,val);
			break;
		case 'l':
			L.setData(data.n,index,val);
			break;
		case 'g':
			G.setData(data.n,index,val);
			break;
		case 'T':
			disT.setData(data.n,index,val);
			break;
		case 'X':
			disX.setData(data.n,index,val);
			break;
		case 'U':
			disU.setData(data.n,index,val);
			break;
	}
}

double DataStorage::getData(const Selector &data, int index) const {

	switch (data.c) {
		case 't':
			return T.getData(data.n,index);
		case 'x':
			return X.getData(data.n,index);
		case 'u':
			return U.getData(data.n,index);
		case 'l':
			return L.getData(data.n,index);
		case 'g':
			return G.getData(data.n,index);
		case 'T':
			return disT.getData(data.n,index);
		case 'X':
			return disX.getData(data.n,index);
		case 'U':
			return disU.getData(data.n,index);
	}
	return 0.0;
}


void DataStorage::freeze() {
	T.freeze();
	X.freeze();
	U.freeze();
	G.freeze();
	L.freeze();
	disT.freeze();
	disX.freeze();
	disU.freeze();
}

void DataStorage::unFreeze() {
	T.unFreeze();
	X.unFreeze();
	U.unFreeze();
	G.unFreeze();
	L.unFreeze();
	disT.unFreeze();
	disX.unFreeze();
	disU.unFreeze();
}

void DataStorage::copy(DataStorage &ds) {
	T.copy(ds.T);
	X.copy(ds.X);
	U.copy(ds.U);
	G.copy(ds.G);
	L.copy(ds.L);
	disT.copy(ds.disT);
	disX.copy(ds.disX);
	disU.copy(ds.disU);
}

double DataStorage::getTF() const {
	return tfgesamt;
}

double DataStorage::getT0() const {
	return t0gesamt;
}

int DataStorage::getLength() const {
	return length;
}

int DataStorage::getDisLengthX() const {
	return disLengthX;
}

int DataStorage::getDisLengthU() const {
	return disLengthU;
}

}
