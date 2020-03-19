//
// C++ Interface: plot
//
// Description:
//
//
// Author: Matthias Knauer <knauer@math.uni-bremen.de>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef plot_h
#define plot_h

#include "../base/point.h"
#include <string>
#include <vector>
#include <cstring>
#include <iostream>
#include "xopt_data.h"
#include "../glbase/smoothmovement.h"
#include "SDL2/SDL.h"
#include "functions.h"

#include "../core/TWparameter.h"

enum cols {Black, White, Red, Grey, Green, Cyan, Rose, Hue=10, Hue2=110};

extern int onlyhigh;

class EpsWriter;
class Viewer;
class SDLFrame;

class BasePlot {
public:
	BasePlot(int ind);
	virtual ~BasePlot();
	virtual void Draw(DataStorage &ds, DataStorage &dstop, int ii);
	virtual void DrawMore(SDLFrame *viewer, double *x, double t) {}


	virtual void DrawText() const;
	virtual void DrawText00(double CurTime, double PlTime) const {}
	virtual void drawData(DataStorage &ds, DataStorage &dstop, int ii) {}
	virtual void SetMaxTime(DataStorage::TimeMode_e timemode, double time) {}
	virtual void SetMaxTime2(double time) {}
	/*	virtual int GetType() const {
			return 0;
		}
	*/
	virtual void SetUserControl(FunktionenzeigerC c) {
		//ctrldata=1;
		//FuncC=c;
	}
	void newHighLow(DataStorage &ds);
	virtual bool MouseInput(DataStorage &ds, int button, const  Point<int> &p) {
		return false;
	}
	void SetColor(cols c) const;
	virtual bool SpecialKeys(const Uint8 *keys) {
		return false;
	}
	virtual void TimerFunc(double t) {}
	virtual void Info() {}

	virtual void DrawCompareCurve(DataStorage &ds) const;
	virtual bool contains(const Point<int> &p) ;
	int containsArea(const Point<int> &p) ;
	void calcYscale();
	bool IsIcon() const;
	static int minwidth,minheight;

	//void AddCompareCurve(double *cmp, int cmpstep);
	void AddCompareCurve2(double *cmptime, double *cmp, int cmpstep, int n);

	void AddDynCompareCurve(double *time, double *cmp, int *cmpstep, int *n);
	void DrawDynCompareCurve(DataStorage &ds) const;

	/*void SetSubname(char* s) {
		strcpy(epssubname,s);
		//epssubname[strl]=0;
	}*/
	void SetEpsTitle(const std::string &buf) {
		epstitle=buf;
	}
	void SetScaleData(double sd) {
		scaledata=sd;
	}
	void ControlData(int i) {
		ctrldata = i;
	}

	void DrawString(int x, int y, const std::string &s) const ;
	void DrawLine(double x1, double y1, double x2, double y2, int mode);
	void DrawPolyline(double *x, double *y, int n, int mode);


	struct Geometry {
		Geometry(int w=0, int h=0);
		Point<int> pos;
		int width;
		int height;
	};

	struct DynGeometry {
		DynGeometry(int w, int h);
		void SetNextGeometry(int currenttime, const Geometry &g, int time);
		bool Update(int &change, int currenttime);
		int lasttime;
		int nexttime;
		Geometry last;
		Geometry next;
		Geometry cur;
	}
	geom;

public:
	void MoveTo(int currenttime, const Geometry &g, int time) {
		geom.SetNextGeometry(currenttime, g,time);
	}
	Point<int> GetPos() const {
		return geom.cur.pos;
	}
	int GetWidth() const {
		return geom.cur.width;
	}
	int GetHeight() const {
		return geom.cur.height;
	}

	void Matlab();
	void Print(DataStorage &ds, DataStorage &dstop);
	
	int redflag;

	/** Left border */
	static int LBorder;
	/** Right border */
	static int RBorder;
	/** Upper border */
	static int TBorder;
	/** Lower border */
	static int BBorder;

	void epsCompareCurve(EpsWriter *epsw) const;
	virtual void epsData(DataStorage &ds, DataStorage &dstop, EpsWriter *epsw)=0;
	virtual void matlab(std::ostream &os) const {};


protected:
	virtual void drawFrameText() const ;
	virtual void drawFrame() const;
	virtual void drawIconFrame() const;
	virtual void drawPlotName() const;

public:
	float MapFloatToX(double d) const;
	virtual float MapToY(double d) const;
	float MapStopTimeToX(DataStorage &ds, DataStorage &dstop, int i) const;
	float MapTimeToX(DataStorage &ds, int i) const;

	void Timer(Viewer *v);
	int extx,exty;

	static int allplotnames;
	static int acsize;
	std::string epstitle;

	virtual double roundStep(double) const;

	std::vector< std::vector<int> > dot_indices; // pro Phase!!!
protected:
	virtual double MapFromX(short x) const;
	double MapFromY(short y) const;

	virtual int rawHighLow(DataStorage &ds)=0;
	int hasControlData() const {
		return ctrldata;
	}
	bool Mirror;
	

public:
	double High,Low,High2, Low2;
protected:
	int index;
	double *comp;
	double *comptime;
	int mouseOver;
	int compn;

	double *dyncomp;
	double *dyncomptime;
	int *dyncompn;
	int *dyncompstep;

	double yscale;
	double scaledata;

	int compstep;
	
	int ctrldata;

	float GetHue() const;

};


class GeneralPlot : public BasePlot {
public:
	GeneralPlot(char c, int d, int ind);
	virtual ~GeneralPlot();
	void drawData(DataStorage &ds, DataStorage &dstop, int ii);
	void SetMaxTime(DataStorage::TimeMode_e timemode, double time);
	
	/*	int GetType() const {
			return 1;
		}
	*/
	int GetDgl() const {
		return data.n;
	}
	virtual bool MouseInput(DataStorage &ds, int button,  const Point<int> &p) ;
	double MapFromX(DataStorage &ds, short x) const;
protected:
	virtual void epsData(DataStorage &ds, DataStorage &dstop, EpsWriter *epsw) ;
	int rawHighLow(DataStorage &ds);

private:
	Selector data;


};


class PhasePlot : public BasePlot {
public:
	PhasePlot(char c1, int d1, char c2, int d2, int ind);
	virtual ~PhasePlot();

protected:
	Selector data1,data2;

	virtual int rawHighLow(DataStorage &ds);
	virtual void drawData(DataStorage &ds, DataStorage &dstop, int ii);
	virtual void epsData(DataStorage &ds, DataStorage &dstop, EpsWriter *epsw) ;

};


class MatrixPlot : public BasePlot {
public:
	MatrixPlot(double *d, int dim1,int dim2,int ind);
	virtual ~MatrixPlot();

protected:
	int D1,D2;
	double *matrix;

	virtual int rawHighLow(DataStorage &ds);
	virtual void drawData(DataStorage &ds, DataStorage &dstop, int ii);
	virtual void epsData(DataStorage &ds, DataStorage &dstop, EpsWriter *epsw);

};


#include "worhp/C_cs.h"
#include <../glbase/globject.h>

class SparsePlot : public BasePlot {
public:
	SparsePlot(WorhpMatrix *d);
	virtual ~SparsePlot();
	WorhpMatrix *matrix;
protected:


	virtual int rawHighLow(DataStorage &ds);
	virtual void drawData(DataStorage &ds, DataStorage &dstop, int ii);
	virtual void epsData(DataStorage &ds, DataStorage &dstop, EpsWriter *epsw);
	virtual void matlab(std::ostream &os) const;

	virtual double roundStep(double) const;
	virtual bool contains(const Point<int> &p) ;
	float MapToY(double d) const;
	std::string infostring;
	virtual void drawPlotName() const;
	Point<int>mouse;
};


class DataPlot : public BasePlot {
public:
	DataPlot(std::string &s, Funktionenzeiger2 func, int par, int ind);
	DataPlot(std::string &s, Funktionenzeiger2 func, int *par, int ind);

	virtual ~DataPlot();
	virtual void SetMaxTime(DataStorage::TimeMode_e timemode, double time);

protected:
	Funktionenzeiger2 Func;
	std::string title;
	//	int param;
	virtual int rawHighLow(DataStorage &ds);
	//void drawData();
	void drawData(DataStorage &ds,DataStorage &dstop, int ii);
	virtual void epsData(DataStorage &ds, DataStorage &dstop, EpsWriter *epsw) ;

	std::vector<int> indices;

};

/**
* Klasse zum Plotten des Fehlers und der Schrittweite
* @author Matthias Rick
*/
class gitterPlot : public BasePlot {
public:
	gitterPlot(std::string &s, int i, const std::vector<std::vector<double> > &zeit, const std::vector<std::vector<double> > &fehler);
	~gitterPlot();

protected:
	std::string title;
	/** Stuetzstellen */
	const std::vector<std::vector<double> > &zeit;
	/** Fehler */
	const std::vector<std::vector<double> > &fehler;

	std::string ticktext2(double step, double i) const;
	void drawFrameText() const;
	virtual int rawHighLow(DataStorage &ds);
	void drawData(DataStorage &ds, DataStorage &dstop, int ii);
	virtual std::string getLabel() const;
	virtual void epsData(DataStorage &ds, DataStorage &dstop, EpsWriter *epsw);
};

/**
* Klasse zum Plotten der Stuetzstellen
* @author Matthias Rick
*/
class punktPlot : public BasePlot {
public:
	/** Erstellt einen Plot der Stuetzstellen
	* @param s Titel
	* @param i
	* @param zeit Diskretisierung
	* @param twdiscretization Diskretisierungstyp
	*/
	punktPlot(std::string &s, int i, const std::vector<std::vector<double> > &zeit, const TWdiscretization *twdiscretization);
	~punktPlot();

protected:
	/** Titel */
	std::string title;
	/** Stuetzstellen */
	const std::vector<std::vector<double> > &zeit;
	/** Diskretisierungstyp */
	const TWdiscretization *twdiscretization;

	virtual int rawHighLow(DataStorage &ds);
	void drawData(DataStorage &ds, DataStorage &dstop, int ii);
	virtual std::string getLabel() const;
	virtual void epsData(DataStorage &ds, DataStorage &dstop, EpsWriter *epsw);
};

/**
* Klasse zum Plotten der Beschraenkungen
* @author Matthias Rick
*/
class lambdaPlot : public BasePlot {
public:
	/** Erstellt einen Plot der Beschraenkungen (Lambda)
	* @param s Titel
	* @param zeit Diskretisierung
	* @param lambda Multiplikatoren
	* @param n_dis Anzahl Stuetzstellen
	* @param n_ctrl Anzahl Steuerungen
	* @param n_ode Anzahl Zustaende
	* @param twdiscretization Diskretisierungstyp
	*/
	lambdaPlot(std::string &s, const double *zeit, const double *lambda, int n_dis, int n_ctrl, int n_ode, const TWdiscretization *twdiscretization);
	~lambdaPlot();

protected:
	/** Titel */
	std::string title;
	/** Stuetzstellen */
	const double *zeit;
	/** Multiplikatoren */
	const double *lambda;
	/** Anzahl Stutzstellen */
	int n_dis;
	/** Anzahl Steuerungen */
	int n_ctrl;
	/** Anzahl Zustaende */
	int n_ode;
	/** Diskretisierungstyp */
	const TWdiscretization *twdiscretization;

	int rawHighLow(DataStorage &ds);
	virtual void drawData(DataStorage &ds, DataStorage &dstop, int ii);
	std::string getLabel() const;
	void epsData(DataStorage &ds, DataStorage &dstop, EpsWriter *epsw);
};

/**
* Klasse zum Plotten einer bestimmten Beschraenkung
* @author Matthias Rick
*/
class lambdaPlot2 : public lambdaPlot {
	
public:
	/** Erstellt einen Plot der Beschraenkungen (Lambda) fuer eine bestimme Komponente
	* @param s Titel
	* @param n Komponente
	* @param type Zustand=0 oder Steuerung=1
	* @param zeit Diskretisierung
	* @param lambda Multiplikatoren
	* @param n_dis Anzahl Stuetzstellen
	* @param n_ctrl Anzahl Steuerungen
	* @param n_ode Anzahl Zustaende
	* @param twdiscretization Diskretisierungstyp
	*/
	lambdaPlot2(std::string &s, int n, int type, const double *zeit, const double *lambda, int n_dis, int n_ctrl, int n_ode, const TWdiscretization *twdiscretization);
	~lambdaPlot2();
	
protected:
	/** Komponente (zB Zustand 3) */
	int komp;
	/** Typ (Zustand=0 oder Steuerung=1) */
	int type;
	/** Vector zum zwischenspeichern der Plot-Daten */
	std::vector<double> plotData;
	
	void drawData(DataStorage &ds, DataStorage &dstop, int ii);
};

/**
* Klasse zum Plotten der Adjungierten (Mu)
* @author Matthias Rick
*/
class adjPlot : public BasePlot {

public:
	/** Erstellt einen Plot der Adjungierten (Mu)
	* @param s Titel
	* @param zeit Diskretisierung
	* @param mu Adjungierten
	* @param n_dis Anzahl Stuetzstellen
	* @param n_ctrl Anzahl Steuerungen
	* @param n_ode Anzahl Zustaende
	* @param n_con 
	* @param twdiscretization Diskretisierungstyp
	*/
	adjPlot(std::string &s, const double *zeit, const double *mu, int n_dis, int n_ctrl, int n_ode, int n_con, const TWdiscretization *twdiscretization);
	~adjPlot();

protected:
	/** Titel */
	std::string title;
	/** Stuetzstellen */
	const double *zeit;
	/** Adjungierten */
	const double *mu;
	/** Anzhal Stuetzstellen */
	int n_dis;
	/** Anzahl Steuerungen */
	int n_ctrl;
	/** Anzahl Zustaende */
	int n_ode;
	/** Anzahl aller Beschraenkungen */
	int n_con;
	/** Diskretisierungstyp */
	const TWdiscretization *twdiscretization;

	virtual void drawData(DataStorage &ds, DataStorage &dstop, int ii);
	std::string getLabel() const;
	void epsData(DataStorage &ds, DataStorage &dstop, EpsWriter *epsw);
	int rawHighLow(DataStorage &ds);
};

/**
* Klasse zum Plotten einer bestimmten Adjungierten
* @author Matthias Rick
*/
class adjPlot2 : public adjPlot {
public:
	/** Erstellt einen Plot einer Adjungierten (Mu)
	* @param s Titel
	* @param n Komponente
	* @param zeit Diskretisierung
	* @param mu Adjungierten
	* @param n_dis Anzahl Stuetzstellen
	* @param n_ctrl Anzahl Steuerungen
	* @param n_ode Anzahl Zustaende
	* @param n_con 
	* @param twdiscretization Diskretisierungstyp
	*/
	adjPlot2(std::string &s, int n, const double *zeit, const double *mu, int n_dis, int n_ctrl, int n_ode, int n_con, const TWdiscretization *twdiscretization);
	~adjPlot2();
	
protected:
	/** Komponente */
	int komp;
	/** Vector zum zwischenspeichern der Plot-Daten */
	std::vector<double> plotData;
	
	void drawData(DataStorage &ds, DataStorage &dstop, int ii);
};


class Data2Plot : public BasePlot {
public:
	Data2Plot(std::string& s, Funktionenzeiger2 func, int par1, int par2, int ind);

	virtual ~Data2Plot();

protected:
	Funktionenzeiger2 Func;
	std::string title;
	//	int param;
	virtual int rawHighLow(DataStorage &ds);
	void drawData(DataStorage &ds, DataStorage &dstop, int ii);
	virtual void epsData(DataStorage &ds, DataStorage &dstop, EpsWriter *epsw) ;

	std::vector<int> indices;

};

class UserPlot : public BasePlot {
public:
	UserPlot(Funktionenzeiger func,  int ind);
	UserPlot(FunktionenzeigerI func, int ind, int index_);
	UserPlot(FunktionenzeigerU func, int ind);
	virtual ~UserPlot();
	/*	int GetType() const {
			return 2;
		}
	*/
	void SetMaxTime2(double time);


	void SetUserControl(FunktionenzeigerC c) {
		ctrldata=1;
		FuncC=c;
	}

	virtual bool MouseInput(DataStorage &ds, int button, const Point<int> &p);

	void CallUserControl(DataStorage &ds, int button, double x, double y);

protected:
	Funktionenzeiger Func;
	FunktionenzeigerI FuncI;
	FunktionenzeigerC FuncC;
	FunktionenzeigerU FuncU;
	double maxtime;
	int lasttime;

	double userBorder[4];

	virtual int rawHighLow(DataStorage &ds);
	virtual void drawData(DataStorage &ds, DataStorage &dstop, int ii);
	
	int iindex;
	virtual void epsData(DataStorage &ds, DataStorage &dstop, EpsWriter *epsw) ;

};

class TabularPlot : public BasePlot {
public:
	TabularPlot(int n, int m, double *val, char *heads, int ind);

	virtual ~TabularPlot();
	//	virtual void SetMaxTime(double time);

protected:
	int nn,mm;
	double *values;
	char *headers;

	std::string title;
	//	int param;
	virtual int rawHighLow(DataStorage &ds);
	void drawData(DataStorage &ds, DataStorage &dstop, int ii) {}
	
	virtual void epsData(DataStorage &ds, DataStorage &dstop, EpsWriter *epsw) {}

	void drawFrame() const;
	void drawIconFrame() const;
	void drawFrameText() const;
	void drawTextData() const;
	void drawIconTextData() const;
};


//class SDLFrame;
class XMLNode;
class glObject;


class Camera {
public:
	Camera() {}
	void Init(XMLNode *n);
	float x;
	float y;
	float z;
	float rotx;
	float roty;
	float focus;

};



typedef void (*plot3d) (glObject *obj, double *x, double t);


class ThreeDPlot : public BasePlot {
public:
	ThreeDPlot(XMLNode *xml, plot3d f);
	~ThreeDPlot();

protected:
	void Draw(DataStorage &ds, DataStorage &dstop);
	
	void epsData(DataStorage &ds, DataStorage &dstop, EpsWriter *epsw) {}
	void drawBack() const;
	int rawHighLow(DataStorage &ds) {
		return 0;
	}
	void DrawText() const {}
	void DrawText00(double CurTime, double PlTime) const;
	//void displayObjects(double *x, double t);
	void drawFrame() const;

	void DrawMore(SDLFrame *viewer, double *x, double t);
	bool SpecialKeys(const Uint8 *keys);
	std::vector<glObject*> obj;
	SmoothMovement campos[6];
	Camera cam[10];
	void TimerFunc(double t);
	void Info();
	plot3d func3d;


};




#endif // plot_h
