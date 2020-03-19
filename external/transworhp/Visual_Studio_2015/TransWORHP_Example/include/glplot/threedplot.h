#pragma once

#include "baseplot.h"

#include "../glbase/smoothmovement.h"

class XMLNode;

namespace tw {

class glObject;

using plot3d = void (*) (glObject *obj, double *x, double t);

struct Camera {
	void Init(XMLNode *n);
	float x;
	float y;
	float z;
	float rotx;
	float roty;
	float focus;
};

class ThreeDPlot : public BasePlot {
public:
	ThreeDPlot(XMLNode *xml, plot3d f);
	~ThreeDPlot();

protected:
	void Draw(DataStorage &ds, DataStorage &dstop);

	void epsData(DataStorage &ds, DataStorage &dstop, EpsWriter *epsw) override;
	void drawBack() const;
	int rawHighLow(DataStorage &ds) override;
	void drawText() const override;
	void DrawText00(double CurTime, double PlTime) const override;
	//void displayObjects(double *x, double t);
	void drawFrame() const override;

	void DrawMore(SDLFrame *viewer, double *x, double t) override;
	bool SpecialKeys(const Uint8 *keys) override;

	std::vector<glObject*> obj;
	SmoothMovement campos[6];
	Camera cam[10];
	void TimerFunc(double t) override;
	void Info() override;
	plot3d func3d;
};

}
