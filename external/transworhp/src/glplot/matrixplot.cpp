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

MatrixPlot::MatrixPlot(double *d, int dim1,int dim2,int ind)
		: BasePlot(ind), D1(dim1), D2(dim2) {

				Low = 0;
				High=dim1;
				
	Low2 = 0;
	High2 = dim2;

	matrix = d;

	/*Dgl = d1;
	Dgl2 = d2;*/
}

MatrixPlot::~MatrixPlot() {}

void MatrixPlot::drawData(DataStorage &/*ds*/, DataStorage &/*dstop*/, int /*ii*/) {
	int curx,cury;//,lastx,lasty;

	glBegin(GL_POINTS);

	for (int i=0;i<D1;i++) {
		for (int j=0;j<D2;j++) {

			cury = (int)MapToY(j);
			curx = (int)MapFloatToX(i);

			double dd = matrix[i+D2*j];

			if (dd>0.001) {
				SetColor(TWcolor::Red);
				glVertex3f((float)curx, (float)cury,0.5f);

			}
			if (dd<-0.001) {
				SetColor(TWcolor::Green);
				glVertex3f((float)curx, (float)cury,0.5);
			}

			/*lastx=curx;
			lasty=cury;*/
		}
	}
	glEnd();
	/*
		if (Dgl>0)
			SetForeground(red);
		else
			SetForeground(white);
	*/
	

}





int MatrixPlot::rawHighLow(DataStorage &/*ds*/) {
	/*double v;
	for (int i=0;i<gl->ds.GetLength();i++) {
		v = gl->ds.GetData(Dgl,i);
		if (v < Low || i==0)
			Low = v;
		if (v > High || i==0)
			High = v;
	}
	for (int i=0;i<gl->ds.GetLength();i++) {
		v = gl->ds.GetData(Dgl2,i);
		if (v < Low2 || i==0)
			Low2 = v;
		if (v > High2 || i==0)
			High2 = v;
	}*/
	return 1;

}


void MatrixPlot::epsData(DataStorage &/*ds*/, DataStorage &/*dstop*/, EpsWriter */*epsw*/) {

/*
	epsw->SetLineColor(Red);// "col2";

	for (int i=0;i<gl->ds.GetLength();i++) {
		epsw->AddPoint(epsw->epsMapFloatToX(gl->ds.GetData(Dgl2,i)),
		               epsw->epsMapToY(gl->ds.GetData(Dgl,i)));
	}

	epsw->Line();
	epsw->ClearPoints();


	if (gl->dstop.GetLength()) {

		epsw->SetLineColor(Green);// "col2";

		for (int i=0;i<gl->dstop.GetLength();i++) {
			epsw->AddPoint(epsw->epsMapFloatToX(gl->dstop.GetData(Dgl2,i)),
			               epsw->epsMapToY(gl->dstop.GetData(Dgl,i)));
		}
		epsw->Line();
		epsw->ClearPoints();
	}*/
}

}
