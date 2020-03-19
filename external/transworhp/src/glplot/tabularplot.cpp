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
#include "conversion.h"
#include "../toolbase/tool.h"

#include <iostream>
#include <iomanip>

#ifdef WIN32
#define fabs(A) ((A)>0?(A):-(A))
#else
#define min(A,B) ((A)<(B)?(A):(B))
#define max(A,B) ((A)>(B)?(A):(B))
#endif

using namespace std;

namespace tw {

TabularPlot::TabularPlot(int n, int m, double *val, char *heads, int ind)
		: BasePlot(ind), nn(n), mm(m), headers(heads) {

			Low=0;
			High=m+1;
			
	High2 = n;
	Low2 = 0;
	values = val;

	if (n>1) extx=2;
}

TabularPlot::~TabularPlot() {}

void TabularPlot::drawTextData() const {


	glColor4f(.8,0,0,1.);
	float yzero=MapToY(mm + 1-.5);

	for (int j=0;j<nn;j++) {


		std::string buf(&headers[j*100],100);
		int ww = Tool::font->StringLength(buf.c_str());

		float xzero=MapFloatToX(j+1-.5);
		int pos = xzero-ww/2;
		DrawString(pos,yzero-5,buf);


	}

	glColor4f(0,0,0,.7);

	for (int i=0;i<mm;i++) {
		float yzero=MapToY(mm - i-.5);

		for (int j=0;j<nn;j++) {
			double a = values[j+i*nn];
			std::string buf = ToString(a);
			int ww = Tool::font->StringLength(buf.c_str());

			float xzero=MapFloatToX(j+1-.5);
			int pos = xzero-ww/2;
			DrawString(pos,yzero-5,buf);


		}
	}
}



void TabularPlot::drawIconTextData() const {

	int j=0;

	glColor4f(0,0,0,.7);

	for (int i=0;i<mm;i++) {
		float yzero=MapToY(mm - i);

		double a = values[j+i*nn];
		std::string buf = ToString(a);
		int ww = Tool::font->StringLength(buf.c_str());

		float xzero=MapFloatToX(nn/2);
		int pos = xzero-ww/2;
		DrawString(pos,yzero-5,buf);

	}
}




int TabularPlot::rawHighLow(DataStorage &/*ds*/) {
	return 0;
}

void TabularPlot::drawData(DataStorage&, DataStorage&, int) {}
	
void TabularPlot::epsData(DataStorage&, DataStorage&, EpsWriter*) {}


/*
void TabularPlot::epsData(EpsWriter *epsw) {
 
	std::vector<int>::iterator it = indices.begin();
	for (;it!=indices.end();it++) {
 
		epsw->SetLineColor(Red);// "col2";
		for (int i=0;i<gl->ds.GetLength();i++) {
			epsw->AddPoint(epsw->epsMapTimeToX(i),
			               epsw->epsMapToY(gl->ds.Call(Func,i,*it)));
		}
		epsw->Line();
		epsw->ClearPoints();
 
 
		if (gl->dstop.GetLength()) {
			epsw->SetLineColor(Green);// "col3";
 
			for (int i=0;i<gl->dstop.GetLength();i++) {
				epsw->AddPoint(epsw->epsMapStopTimeToX(i),
				               epsw->epsMapToY(gl->dstop.Call(Func,i,*it)));
			}
			epsw->Line();
			epsw->ClearPoints();
		}
	}
}
*/



void TabularPlot::drawFrame() const {

	int width = GetWidth();
	int height = GetHeight();

	glBegin(GL_QUADS);
	float a = GetHue()*.2+.3; // [.3 ; .5]

	int z1 = (width-RBorder-LBorder)*.2;
	int z2 = (height-TBorder-BBorder)*.2;

	int z = min(z1,z2);
	//z = min(z1,10);

	glColor4f(1,1,1,a);
	glVertex3f(LBorder,TBorder,0);
	glVertex3f(width-RBorder,TBorder,0);
	glColor4f(1,1,1,a+.4);
	glVertex3f(width-RBorder,height-BBorder,0);
	glVertex3f(LBorder,height-BBorder,0);
	glColor4f(1,1,1,a);
	glVertex3f(LBorder,TBorder,0);
	glEnd();



	/*double step = 25*(High2-Low2)/(width-RBorder-LBorder);
	step = roundStep(step);
	*/

	float col[3][4] = {{0.f,0.f,0.f,.4f},{0.f,0.f,0.f,.7f},{0.f,0.f,0.f,.1f}};

	for (double i=0;i<=nn;i++) {
		float xzero=MapFloatToX(i);

		glBegin(GL_LINE_STRIP);
		if (xzero>LBorder+2 && xzero<width-RBorder-2) {
			glColor4fv(col[i>0?0:1]);
			glVertex3f(xzero, TBorder,0.1);
			glColor4fv(col[i>0?2:0]);
			glVertex3f(xzero, TBorder+z,0.1);
			glVertex3f(xzero, height-BBorder-z,0.1);
			glColor4fv(col[i>0?0:1]);
			glVertex3f(xzero,height-BBorder,0.1);
		}
		glEnd();

	}


	/*
		glColor4f(0,0,.9,1);
		step = 15*(High-Low)/(height-TBorder-BBorder);
		step = roundStep(step);
	 
		double l = max(0,Low);
		l = (int)(l/step)*step;*/

	for (double i=0;i<=mm;i++) {
		float yzero=MapToY(i);

		glBegin(GL_LINE_STRIP);

		if (yzero>TBorder+2 && yzero<height-BBorder-2) {
			glColor4fv(col[i>0?0:1]);
			glVertex3f(LBorder,yzero,0.1);
			glColor4fv(col[i>0?2:0]);
			glVertex3f(LBorder+z,yzero,0.1);
			glVertex3f(width-RBorder-z,yzero,0.1);
			glColor4fv(col[i>0?0:1]);
			glVertex3f(width-RBorder,yzero,0.1);
		}

		glEnd();

	}




	glBegin(GL_LINE_STRIP);
	glColor4f(0,0,0,.5);
	glVertex3f(LBorder,TBorder,0.1);
	glVertex3f(width-RBorder,TBorder,0.1);
	glVertex3f(width-RBorder,height-BBorder,0.1);
	glVertex3f(LBorder,height-BBorder,0.1);
	glVertex3f(LBorder,TBorder,0.1);
	glEnd();


	glBegin(GL_LINE_STRIP);
	glColor4f(0,0,0,.1);
	glVertex3f(0,0,0.1);
	glVertex3f(0,height,0.1);
	glVertex3f(width,height,0.1);
	glVertex3f(width,0,0.1);
	glVertex3f(0,0,0.1);
	glEnd();

	drawTextData();
}




void TabularPlot::drawIconFrame() const {

	int width = GetWidth();
	int height = GetHeight();

	glBegin(GL_QUADS);
	float a = GetHue()*.2+.3; // [.3 ; .5]


//	int z1 = (width-0-0)*.2;
//	int z2 = (height-0-0)*.2;

	//int z = min(z1,z2);
	//z = min(z1,10);

	glColor4f(1,1,1,a);
	glVertex3f(0,0,0);
	glVertex3f(width,0,0);
	glColor4f(1,1,1,a+.4);
	glVertex3f(width,height,0);
	glVertex3f(0,height,0);
	glColor4f(1,1,1,a);
	glVertex3f(0,0,0);
	glEnd();


	
	glBegin(GL_LINE_STRIP);
	glColor4f(0,0,0,.5);
	glVertex3f(0,0,0.1);
	glVertex3f(width-0,0,0.1);
	glVertex3f(width-0,height-0,0.1);
	glVertex3f(0,height-0,0.1);
	glVertex3f(0,0,0.1);
	glEnd();

	drawIconTextData();
}


void TabularPlot::drawFrameText() const {

	int width = GetWidth();
	//int height = GetHeight();
	double step = 25*(High2-Low2)/(width-RBorder-LBorder);
	step = roundStep(step);

	for (int i=0;i<mm;i++) {
		float yzero=MapToY(mm - i-.5);

		std::string buf(&headers[(i+nn)*100],100);


		int ww = Tool::font->StringLength(buf.c_str());
		int pos = LBorder-ww-2;

		DrawString(pos,yzero-5,buf);
	}
}

}
