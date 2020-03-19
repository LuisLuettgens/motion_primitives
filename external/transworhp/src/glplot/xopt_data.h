#pragma once

#include "functions.h"

namespace tw {

struct Selector {

	Selector(char c_, int n_) : c(c_),n(n_) {}

	char c; // 't', 'x', 'u', 'l', 'g'
	int n;  // 0..
};



/** @defgroup data Data Access and Storage
 *  ???
 * @{
 */
class DataStorage;

/** Treat FORTRAN vector as 2-dimensional vector */
class DataVector {
	friend class DataStorage;
public:
	/** Constructor */
	DataVector();
	/** Destructor */
	~DataVector() {}

	void setVector(int n, int s, int off, double *value);
	double getData(int dgl, int index) const;
	void setData(int dgl, int index, double v);
	void freeze();
	void unFreeze();
	void copy(DataVector &other);
	int getSize() const;

//private:

	/** Anzahl der diskreten Daten */
	int ndis;

	/** Werte pro diskretem Punkt */
	int size;

	/** SpeicherOffset zu naechsten Datenpunkt */
	int offset;

	/** Daten */
	double *D;
};

/** Contains vectors for time, x, u, g */
class DataStorage {
public:
	DataStorage();
	~DataStorage() {}

	void setData(int len, double *t,
			double *x, int xoffset, int xdata,
			double *u, int uoffset, int udata,
			double *g, int goffset, int gdata,
			double *l, int loffset, int ldata);

	void addDisData(int disLenX, int disLenU,
			
			double *dis_t,
			double *dis_x, double *dis_u);
	
	void addDisData(int disLenX, int disLenU,
			int xoffset, int uoffset,
			double *dis_t,
			double *dis_x, double *dis_u);

	void setData(const Selector &data,int index, double val);
	double getData(const Selector &data,int index) const;

	void freeze();
	void unFreeze();
	void copy(DataStorage &ds);


	double getTF() const;
	double getT0() const;
	int getLength() const;
	int getDisLengthX() const;
	int getDisLengthU() const;

	void call(Funktionenzeiger Func, int index, int f, double *bord) {
		(*Func)(index,f,length,X.size,U.size,T.D,X.D,U.D,bord);
	}
	void call(FunktionenzeigerI Func, int index, int f, double *bord,int iindex) {
		(*Func)(index,f,length,X.size,U.size,T.D,X.D,U.D,bord,iindex);
	}
	double call(Funktionenzeiger2 Func, int i, int param) {
		return (*Func)(length,X.size,U.size,T.D,X.D,U.D,i,param);
	}
	void call(FunktionenzeigerC Func, int button, double xpos, double ypos) {
		(*Func)(button,xpos,ypos,length,X.size,U.size,T.D,X.D,U.D);
	}
	void call(FunktionenzeigerU Func, int index, int f, double *bord) {
		(*Func)(index,f,length,X.size,U.size,UNKNOWN.size,T.D,X.D,U.D,UNKNOWN.D,bord);
	}

	enum TimeMode_e {index_as_time, dgl_as_time, fixed_float_as_time, dgl2_as_time, time_is_time,time_as_func};

	static TimeMode_e timemode;
	static int timedgl;
	static double floattime;
	static double floatstart;

	static int nsplitdis;
	static int nsplitlen;
	static int nsplitdis2;
	static int *nsplitindex;
	static double MaxTF;
	static double MinT0;

	/** Steuerung */
	DataVector U;
	/** Zustand */
	DataVector X;
	/** Zeit */
	DataVector T;
	/** Beschraenkungen(?) nicht genutzt */
	DataVector G;
	/** Integral */
	DataVector L;
	/** ? */
	DataVector UNKNOWN;

	/** diskrete Zeit (fuer pmTW) */
	DataVector disT;
	/** diskreter Zustand (fuer pmTW) */
	DataVector disX;
	/** diskrete Steuerung (fuer pmTW) */
	DataVector disU;

	double start;

private:
	double tfgesamt;
	double t0gesamt;
	int length;
	int disLengthX, disLengthU;
	//int nunbe;
};


double CallTimeFunc(int &i,double *X, double *UNBE, double *T, int &len, int &ndgl, int &nunbe);

//extern DataStorage ds;
//extern DataStorage dstop;

/** @} */

}
