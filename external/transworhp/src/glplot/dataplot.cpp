#include "plot.h"

#include "xopt_eps.h"

#include <GL/gl.h>

#include <iostream>
#include <iomanip>

using namespace std;

namespace tw {

DataPlot::DataPlot(std::string &s, Funktionenzeiger2 func,
                   int par, int ind)
		: BasePlot(ind), Func(func), title(s) {

	indices.push_back(par);
}


DataPlot::DataPlot(std::string &s, Funktionenzeiger2 func,
                   int* par, int ind)
		: BasePlot(ind), Func(func), title(s) {

	int *i = par;
	while (*i) {
		indices.push_back(*i);
		i++;
	}
}

Data2Plot::Data2Plot(std::string &s, Funktionenzeiger2 func,
                   int par1, int par2, int ind)
		: BasePlot(ind), Func(func), title(s) {

	indices.push_back(par1);
	indices.push_back(par2);

	Low2=-10;
	High2=10;

	Low=-1;
	High=1;

}


DataPlot::~DataPlot() {}

Data2Plot::~Data2Plot() {}


void DataPlot::drawData(DataStorage &ds,DataStorage &dstop, int /*ii*/) {
	int curx,cury;//,lastx,lasty;

	DrawCompareCurve(ds);
	DrawDynCompareCurve(ds);

	SetColor(TWcolor::Black);

	glLineWidth(1.5);

	auto it = indices.begin();
	for (;it!=indices.end();it++) {

		//		SetForeground(green);

		glBegin(GL_LINE_STRIP);

		for (int i=0;i<ds.getLength();i++) {
			cury = (int)MapToY(ds.call(Func,i,*it));
			curx = (int)MapTimeToX(ds,i);

			glVertex3f((float)curx, (float)cury,0.5);

	/*		lastx=curx;
			lasty=cury;
	*/	}
		glEnd();
		/*		if (Dgl>0)
					SetForeground(red);
				else
					SetForeground(white);
		*/
		glBegin(GL_LINE_STRIP);
		for (int i=0;i<dstop.getLength();i++) {
			cury = (int)MapToY(dstop.call(Func,i,*it));
			curx = (int)MapStopTimeToX(ds,dstop,i);

			glVertex3f((float)curx, (float)cury,0.5);

	/*		lastx=curx;
			lasty=cury;
	*/	}
		glEnd();
	}
		glLineWidth(1.);

}


void Data2Plot::drawData(DataStorage &ds, DataStorage &dstop, int /*ii*/) {
	int curx,cury;//,lastx,lasty;

	DrawCompareCurve(ds);
	DrawDynCompareCurve(ds);

	SetColor(TWcolor::Black);

	glLineWidth(1.5);

	auto it = indices.begin();
	auto it2 = it;
	it2++;
	//for (;it!=indices.end();it++) {

		//		SetForeground(green);

		glBegin(GL_LINE_STRIP);

		for (int i=0;i<ds.getLength();i++) {
			cury = (int)MapToY(ds.call(Func,i,*it));
			curx = (int)MapFloatToX(ds.call(Func,i,*it2));


		glVertex3f((float)curx, (float)cury,0.5);
		/*lastx=curx;
		lasty=cury;*/



	/*		lastx=curx;
			lasty=cury;
	*/	}
		glEnd();
		/*		if (Dgl>0)
					SetForeground(red);
				else
					SetForeground(white);
		*/
		glBegin(GL_LINE_STRIP);
		for (int i=0;i<dstop.getLength();i++) {
			cury = (int)MapToY(dstop.call(Func,i,*it));
			curx = (int)MapStopTimeToX(ds,dstop,(int)dstop.call(Func,i,*it2));

			glVertex3f((float)curx, (float)cury,0.5);

	/*		lastx=curx;
			lasty=cury;
	*/	}
		glEnd();
	//}
		glLineWidth(1.f);

}




int DataPlot::rawHighLow(DataStorage &ds) {

	bool s = false;

	auto it = indices.begin();
	for (;it!=indices.end();it++) {
		for (int i=0;i<ds.getLength();i++) {
			double v = ds.call(Func,i,*it);
			if (!s || v < Low)
				Low = v;
			if (!s || v > High)
				High = v;
			s = true;
		}
	}
	return 0;
}

int Data2Plot::rawHighLow(DataStorage &ds) {


	auto it = indices.begin();
	    for (int i=0;i<ds.getLength();i++) {
			double v = ds.call(Func,i,*it);
			if (v < Low || i==0)
				Low = v;
			if (v > High || i==0)
				High = v;

		}

		it++;

	for (int i=0;i<ds.getLength();i++) {
			double v = ds.call(Func,i,*it);
		if (v < Low2 || i==0)
			Low2 = v;
		if (v > High2 || i==0)
			High2 = v;
	}

	/*Low=-10;
	High=10;*/

	/*Low2=-10;
	High2=10;*/
	return 1;



}


void DataPlot::SetMaxTime(DataStorage::TimeMode_e timemode, double time) {
	High2 = time;
	if (timemode==DataStorage::index_as_time) {
		High2 = time-1;
	}
}



void DataPlot::epsData(DataStorage &ds, DataStorage &dstop, EpsWriter *epsw) {

	for (auto it = indices.begin(); it != indices.end(); ++it) {

		epsw->SetLineColor(TWcolor::Red);// "col2";
		for (int i=0;i<ds.getLength();i++) {
			epsw->AddPoint(epsw->epsMapTimeToX(ds,i),
			               epsw->epsMapToY(ds.call(Func,i,*it)));
		}
		epsw->Line();
		epsw->ClearPoints();


		if (dstop.getLength()) {
			epsw->SetLineColor(TWcolor::Green);// "col3";

			for (int i=0;i<dstop.getLength();i++) {
				epsw->AddPoint(epsw->epsMapStopTimeToX(ds,dstop,i),
				               epsw->epsMapToY(dstop.call(Func,i,*it)));
			}
			epsw->Line();
			epsw->ClearPoints();
		}
	}
}

void Data2Plot::epsData(DataStorage &ds, DataStorage &dstop, EpsWriter *epsw) {

	for (auto it = indices.begin(); it != indices.end(); ++it) {

		epsw->SetLineColor(TWcolor::Red);// "col2";
		for (int i=0;i<ds.getLength();i++) {
			epsw->AddPoint(epsw->epsMapTimeToX(ds,i),
			               epsw->epsMapToY(ds.call(Func,i,*it)));
		}
		epsw->Line();
		epsw->ClearPoints();


		if (dstop.getLength()) {
			epsw->SetLineColor(TWcolor::Green);// "col3";

			for (int i=0;i<dstop.getLength();i++) {
				epsw->AddPoint(epsw->epsMapStopTimeToX(ds,dstop,i),
				               epsw->epsMapToY(dstop.call(Func,i,*it)));
			}
			epsw->Line();
			epsw->ClearPoints();
		}
	}
}

}
