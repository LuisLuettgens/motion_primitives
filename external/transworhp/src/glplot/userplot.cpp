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

int onlyhigh = 0;



/*void userdraw2_ (int &nhandle, int &nframe,
                 int &ndis1, int &ndgl, int &nsteuer, double *t, double *x,
                 double *u, double *bord);
*/
UserPlot::UserPlot(Funktionenzeiger func, int ind)
		: BasePlot(ind), Func(func), FuncI(nullptr), FuncU(nullptr) {
	/*
		FuncI=0;*/
	
	maxtime = -1;
	lasttime = -1;
}


UserPlot::UserPlot(FunktionenzeigerI func, int ind, int index_)
		: BasePlot(ind), Func(nullptr), FuncI(func), FuncU(nullptr), iindex(index_) {

	/*	Func=0;*/

	maxtime = -1;
	lasttime = -1;
}

UserPlot::UserPlot(FunktionenzeigerU func, int ind)
		: BasePlot(ind), Func(nullptr), FuncI(nullptr), FuncU(func) {
	/*
		FuncI=0;*/

	maxtime = -1;
	lasttime = -1;
}

UserPlot::~UserPlot() {}


void UserPlot::CallUserControl(DataStorage &ds, int button, double x, double y) {
	double xx = x*(High2-Low2)+Low2;
	//(*FuncC)(button,xx,y);
	ds.call(FuncC, button, xx, y);

}


void UserPlot::drawData(DataStorage &ds, DataStorage &/*dstop*/, int /*ii*/) {

	DrawCompareCurve(ds);
	DrawDynCompareCurve(ds);
/*
	int f = gl->GetFrame();

	if (maxtime>=0) {

		if (lasttime>-1) {

			f = (gl->Time-lasttime)/10.;

			if (f>maxtime*100.) {
				lasttime = gl->Time;
			}

		} else {
			lasttime = gl->Time;
		}

	}

	/ *if (dpy->GetTime()>0)
		f = (int)(dpy->GetTime()*1000.);
	* /

	glLineWidth(1.5f);
	if (Func) {
		ds.Call(Func,index,f,userBorder);
	} else if (FuncU) {
		ds.Call(FuncU,index,f,userBorder);
	} else {
		ds.Call(FuncI,index,f,userBorder,iindex);
	}
	if (gl->dstop.GetLength()) {
		/ *if (Dgl>0)
			SetForeground(red);
		else
			SetForeground(white);* /
		if (Func) {
			dstop.Call(Func,index,f,userBorder);
		} else if (FuncU) {
			dstop.Call(FuncU,index,f,userBorder);
		} else {
			dstop.Call(FuncI,index,f,userBorder,iindex);
		}
	}
*/
	glLineWidth(1.f);
}




int UserPlot::rawHighLow(DataStorage &ds) {

	userBorder[0]=1000.;
	userBorder[2]=1000.;
	userBorder[1]=-1000.;
	userBorder[3]=-1000.;

	onlyhigh = 1;

	int f = 0;//gl->GetFrame();
	if (Func) {
		ds.call(Func,index,f,userBorder);
	} else if (FuncU) {
		ds.call(FuncU,index,f,userBorder);
	} else {
		ds.call(FuncI,index,f,userBorder,iindex);
	}

	onlyhigh = 0;
	Low2 = userBorder[0];
	High2 = userBorder[1];
	Low = userBorder[2];
	High = userBorder[3];
	return 1;
}


bool UserPlot::MouseInput(DataStorage &ds, int button, const Point<int> &p) {

	if (hasControlData()) {
		double xx = MapFromX(p.x-GetPos().x);
		double yy = MapFromY(p.y-GetPos().y);
		CallUserControl(ds, button,xx,yy);
		return true;
	}
	return false;
}



void UserPlot::epsData(DataStorage &ds, DataStorage &dstop, EpsWriter *epsw) {

	globEpswriter = epsw;

	epsw->SetLineColor(TWcolor::Red);//="col2";
	int f = 0;//dpy->GetFrame();

	/*if (dpy->GetTime()>0)
		f = (int)(dpy->GetTime()*1000.);
	*/
	if (Func)
		ds.call(Func,index,f,userBorder);
	else if (FuncU)
		ds.call(FuncU,index,f,userBorder);
	else
	    ds.call(FuncI,index,f,userBorder,iindex);

	if (dstop.getLength()) {
		epsw->SetLineColor(TWcolor::Grey);//epscol="col3";
		if (Func)
			dstop.call(Func,index,f,userBorder);
		else if (FuncU)
			dstop.call(FuncU,index,f,userBorder);
		else
			dstop.call(FuncI,index,f,userBorder,iindex);
	}

	globEpswriter = nullptr;

}

void UserPlot::SetUserControl(FunktionenzeigerC c) {
	ctrldata=1;
	FuncC=c;
}

void UserPlot::SetMaxTime2(double time) {

	cout << "SETMAXTIME " << time << endl;
	maxtime = time;
}

}
