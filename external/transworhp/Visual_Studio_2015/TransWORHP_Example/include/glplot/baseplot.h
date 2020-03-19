#pragma once

#include "xopt_data.h"

#include "../base/point.h"

#include <SDL2/SDL.h>

#include <string>
#include <vector>

#ifndef DllExport
#ifdef _WIN32
#define DllExport __declspec(dllexport)
#else
#define DllExport
#endif
#endif

namespace tw {

enum class TWcolor {Black, White, Red, Grey, Green, Cyan, Rose, Hue=10, Hue2=110};

extern int onlyhigh;

class EpsWriter;
class Viewer;
class SDLFrame;

class DllExport BasePlot {
public:
	BasePlot(int ind);
	virtual ~BasePlot();

	void Draw(DataStorage &ds, DataStorage &dstop, int ii);
	virtual void DrawMore(SDLFrame *viewer, double *x, double t);

	virtual void drawText() const;
	virtual void DrawText00(double CurTime, double PlTime) const;
	virtual void drawData(DataStorage &ds, DataStorage &dstop, int ii);
	virtual void SetMaxTime(DataStorage::TimeMode_e timemode, double time);
	virtual void SetMaxTime2(double time);
	virtual void SetMinTime(DataStorage::TimeMode_e timemode, double time);

	virtual void SetUserControl(FunktionenzeigerC c);

	void newHighLow(DataStorage &ds);

	virtual bool MouseInput(DataStorage &ds, int button, const  Point<int> &p);
	virtual bool MouseWheelInput(int scale, const  Point<int>&);

	void SetColor(TWcolor c) const;

	virtual bool SpecialKeys(const Uint8 *keys);

	virtual void TimerFunc(double t);
	virtual void Info();

	virtual void DrawCompareCurve(DataStorage &ds) const;
	virtual bool contains(const Point<int> &p);

	int containsArea(const Point<int> &p);
	void calcYscale();
	/* Fenster ist in Auswahlliste */
	bool IsIcon() const;

	//void AddCompareCurve(double *cmp, int cmpstep);
	void AddCompareCurve2(double *cmptime, double *cmp, int cmpstep, int n);

	void AddDynCompareCurve(double *time, double *cmp, int *cmpstep, int *n);
	void DrawDynCompareCurve(DataStorage &ds) const;

	/*void SetSubname(char* s) {
		strcpy(epssubname,s);
		//epssubname[strl]=0;
	}*/
	void SetEpsTitle(std::string buf);
	void SetScaleData(double sd);
	void ControlData(int i);

	void DrawString(int x, int y, const std::string &s) const;

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

	void MoveTo(int currenttime, const Geometry &g, int time);
	Point<int> GetPos() const;
	int GetWidth() const;
	int GetHeight() const;

	void Matlab();
	void Print(DataStorage &ds, DataStorage &dstop);


	void epsCompareCurve(EpsWriter *epsw) const;
	virtual void epsData(DataStorage &ds, DataStorage &dstop, EpsWriter *epsw)=0;
	virtual void matlab(std::ostream &os) const;

	float MapFloatToX(double d) const;
	virtual float MapToY(double d) const;
	float MapStopTimeToX(DataStorage &ds, DataStorage &dstop, int i) const;
	float MapTimeToX(DataStorage &ds, int i) const;

	void Timer(Viewer *v);

	virtual double roundStep(double) const;


protected:
	virtual void drawFrameText() const;
	virtual void drawFrame() const;
	virtual void drawIconFrame() const;
	virtual void drawPlotName() const;

	double MapFromX(short x) const;
	double MapFromY(short y) const;

	virtual int rawHighLow(DataStorage &ds)=0;
	int hasControlData() const;

	float GetHue() const;

public:

	static int minwidth;
	static int minheight;

	/** Left border */
	static int LBorder;
	/** Right border */
	static int RBorder;
	/** Upper border */
	static int TBorder;
	/** Lower border */
	static int BBorder;

	static int allplotnames;
	static int acsize;

	int redflag;

	int extx, exty;

	std::string epstitle;

	/** Markierungen fuer fixe Anfangs- und Endwerte.
	 * pro Phase */
	std::vector<std::vector<int>> dot_indices;

	// max y-Data
	double High;
	// min y-Data
	double Low;
	//max x-Data
	double High2;
	//min x-Data
	double Low2;

protected:

	int index;
	int mouseOver;
	int compn;
	int compstep;
	int ctrldata;

	double *comp;
	double *comptime;

	double *dyncomp;
	double *dyncomptime;

	int *dyncompn;
	int *dyncompstep;

	double yscale;
	double scaledata;

	double zoomFactorX, zoomFactorY;

	bool Mirror;
};

}
