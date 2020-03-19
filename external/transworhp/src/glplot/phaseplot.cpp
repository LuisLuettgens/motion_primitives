//
// C++ Implementation: plot
//
// Description:
//
//
// Author: Matthias Knauer <knauer@math.uni-bremen.de>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifdef WIN32
#include "windows.h"
#endif
#include "plot.h"
//#include "gldisplay.h"
#include <GL/gl.h>

#include "xopt_eps.h"

#include <iostream>
#include <iomanip>


using namespace std;

namespace tw {

PhasePlot::PhasePlot(char c1, int d1, char c2, int d2,  int ind)
		: BasePlot(ind), data1(c1,d1), data2(c2,d2) {

	//Low2 = lo2;
	//High2 = hi2;
	/*Dgl = d1;
	Dgl2 = d2;*/
}

PhasePlot::~PhasePlot() {}

void PhasePlot::drawData(DataStorage &ds, DataStorage &dstop, int /*ii*/) {
	int curx,cury;//,lastx,lasty;

	glBegin(GL_LINE_STRIP);

	for (int i=0;i<ds.getLength();i++) {
		cury = MapToY(ds.getData(data1,i));
		curx = MapFloatToX(ds.getData(data2,i));

		glVertex3f(curx, cury,0.5);
		/*lastx=curx;
		lasty=cury;*/
	}
	glEnd();
	/*
		if (Dgl>0)
			SetForeground(red);
		else
			SetForeground(white);
	*/
	glBegin(GL_LINE_STRIP);

	for (int i=0;i<dstop.getLength();i++) {
		cury = MapToY(dstop.getData(data1,i));
		curx = MapFloatToX(dstop.getData(data2,i));

		glVertex3f(curx, cury,0.5);

		/*lastx=curx;
		lasty=cury;*/
	}
	glEnd();

}



int PhasePlot::rawHighLow(DataStorage &ds) {
	double v;
	for (int i=0;i<ds.getLength();i++) {
		v = ds.getData(data1,i);
		if (v < Low || i==0)
			Low = v;
		if (v > High || i==0)
			High = v;
	}
	for (int i=0;i<ds.getLength();i++) {
		v = ds.getData(data2,i);
		if (v < Low2 || i==0)
			Low2 = v;
		if (v > High2 || i==0)
			High2 = v;
	}
	return 1;

}


void PhasePlot::epsData(DataStorage &ds, DataStorage &dstop, EpsWriter *epsw) {

	epsw->SetLineColor(TWcolor::Red);// "col2";

	for (int i=0;i<ds.getLength();i++) {
		epsw->AddPoint(epsw->epsMapFloatToX(ds.getData(data2,i)),
		               epsw->epsMapToY(ds.getData(data1,i)));
	}

	epsw->Line();
	epsw->ClearPoints();


	if (dstop.getLength()) {

		epsw->SetLineColor(TWcolor::Green);// "col2";

		for (int i=0;i<dstop.getLength();i++) {
			epsw->AddPoint(epsw->epsMapFloatToX(dstop.getData(data2,i)),
			               epsw->epsMapToY(dstop.getData(data1,i)));
		}
		epsw->Line();
		epsw->ClearPoints();
	}
}

}
