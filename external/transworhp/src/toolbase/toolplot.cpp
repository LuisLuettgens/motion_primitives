//
// C++ Implementation: toolbox
//
// Description:
//
//
// Author: Matthias Knauer,,, <tulio@visurgis>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifdef WIN32
#include "windows.h"
#endif
#include "toolplot.h"
#include "toolview.h"

#include <iomanip>
#include "GL/gl.h"
#include "conversion.h"
#include "toollist.h"
#include "../base/vectortools.h"
#include "../base/status.h"

using namespace std;


#define MAX(A,B) ((A)>(B)?(A):(B))
#define MIN(A,B) ((A)<(B)?(A):(B))
#define ABS(A) ((A)>=0?(A):-(A))


ToolPlot::ToolPlot(const std::string &t, int w, int h, int &items)
		: ToolWindow(t),itemsel(items) {

	LBorder=40;
	RBorder=8;
	TBorder=7+16;
	BBorder=15;


	moving = 0;
	width = 20;
	width = 180+LBorder+RBorder;

	height = 140;
	x = (w - width)
		/2;
	y = (h - height)/2;

	open = 1;
	theight = 18;

	viewit = views.begin();
	currentView = 0;

	AddIcon(ToolIcon::CLOSE);
	AddIcon(ToolIcon::WRITEPNG);
	AddIcon(ToolIcon::WRITETXT);

	resizeable = 3;
	minsizex=100;
	minsizey=100;


}



ToolPlot::~ToolPlot() {
	clearpointer(views);
	ClearPointer(items);

}


void ToolPlot::WriteImage() {
	Message(1000);
}


void ToolPlot::WriteText() {

	if (!currentView)
		return;

	stringstream a;
	a << Tool::filename_prefix + string("_") + currentView->filename + ".dat";

	ofstream of(a.str().c_str());

	of << left;

	stringstream st;
	st << a.str();
	MyStatus("Writing", st.str(), Status::NORMAL);

	currentView->Text(of);

}


inline bool is_round(double value) {
	double v = value-(int)value;
	if (v<0)
		v=-v;
	return (v<1e-8);
}

std::string ticktext(double i) {
	std::stringstream s;
	if (is_round(i)) {
		s << int(i);
	}
	else if (is_round(i*10))
		s << std::setw(1) << i;
	else
		s << std::setw(2) << i;

	return s.str();
}

/** Create Tick text (y) */
std::string ticktext2(double step, double i) {
	std::stringstream s;
	if (step>=1)
		s << std::setw(5) << (int)i;
	else if (step<.1)
		s << std::setiosflags(std::ios::showpoint | std::ios::right |
							  std::ios::fixed)
			<< std::setprecision(2) << std::setw(5) << i;
	else
		s<< setiosflags(std::ios::showpoint | std::ios::right | std::ios::fixed)
			<< std::setprecision(1) << std::setw(5) << i;

	return s.str();
}
double roundStep(double step) {
	if (step<.01)
		step=.01;
	else if (step<.02)
		step=.02;
	else if (step<.05)
		step=.05;
	else if (step<.1)
		step=.1;
	else if (step<.2)
		step=.2;
	else if (step<.5)
		step=.5;
	else if (step<1.)
		step=1.;
	else if (step<2.)
		step=2.;
	else if (step<5.)
		step=5.;
	else if (step<10.)
		step=10.;
	else if (step<20.)
		step=20.;
	else if (step<50.)
		step=50.;
	else if (step<100.)
		step=100.;
	else if (step<200.)
		step=200.;
	else if (step<500.)
		step=500.;
	else if (step<1000.)
		step=1000.;
	else if (step<2000.)
		step=2000.;
	else if (step<5000.)
		step=5000.;
	else if (step<10000.)
		step=10000.;
	else if (step<20000.)
		step=20000.;
	else if (step<50000.)
		step=50000.;
	else if (step<100000.)
		step=100000.;
	else if (step<200000.)
		step=200000.;
	else if (step<500000.)
		step=500000.;
	else if (step<1000000.)
		step=1000000.;
	else if (step<2000000.)
		step=2000000.;
	else /*if (step<5000000.)*/
		step=5000000.;



	return step;
}

void ToolPlot::drawFrame() const {

	// int width = GetWidth();
	int height = this->height - theight;

	glBegin(GL_QUADS);
	float a = .3; // [.3 ; .5]

	// int z1 = (int)((width-RBorder-LBorder)*.2);
	//int z2 = (int)((height-TBorder-BBorder)*.2);

	int z = 5;//min(z1,z2);
	//z = min(z1,10);

	glColor4f(1,1,1,a);
	glVertex3f(LBorder,-height+BBorder,.1);
	glVertex3f(width-RBorder,-height+BBorder,.1);
	glColor4f(1,1,1,a+.4);
	glVertex3f(width-RBorder,-TBorder,.1);
	glVertex3f(LBorder,-TBorder,.1);
	//glColor4f(1,1,1,a);
	//glVertex3f(LBorder,-TBorder,0);
	glEnd();


	double step = getStep();


	float col[3][4] = {{0.f,0.f,0.f,.4f},{0.f,0.f,0.f,.7f},{0.f,0.f,0.f,.1f}};

	for (double i=0;i<=currentView->xmax+.001;i+=step) {
		float xzero=MapFloatToX(i);

		glBegin(GL_LINE_STRIP);
		if (xzero>LBorder+2 && xzero<width-RBorder-2) {
			glColor4fv(col[i>0?0:1]);
			glVertex3f(xzero, -TBorder,0.2);
			glColor4fv(col[i>0?2:0]);
			glVertex3f(xzero, -TBorder-z,0.2);
			glVertex3f(xzero, -height+BBorder+z,0.2);
			glColor4fv(col[i>0?0:1]);
			glVertex3f(xzero,-height+BBorder,0.2);
		}
		glEnd();

	}

	for (double i=-step;i>=currentView->xmin-.001;i-=step) {
		float xzero=MapFloatToX(i);

		glBegin(GL_LINE_STRIP);

		if (xzero>LBorder+2 && xzero<width-RBorder-2) {
			glColor4fv(col[0]);
			glVertex3f(xzero, -TBorder,0.2);
			glColor4fv(col[2]);
			glVertex3f(xzero, -TBorder-z,0.2);
			glVertex3f(xzero, -height+BBorder+z,0.2);
			glColor4fv(col[0]);
			glVertex3f(xzero,-height+BBorder,0.2);
		}
		glEnd();

	}




	glColor4f(0,0,.9,1);
	step = 15*(High-Low)/(height-TBorder-BBorder);
	step = roundStep(step);

	double l = MAX(0,Low);
	l = (int)(l/step)*step;

	for (double i=l;i<=High+.001;i+=step) {
		float yzero=MapToY(i);

		glBegin(GL_LINE_STRIP);

		if (yzero>TBorder+2 && yzero<height-BBorder-2) {
			glColor4fv(col[i>0?0:1]);
			glVertex3f(LBorder,-yzero,0.2);
			glColor4fv(col[i>0?2:0]);
			glVertex3f(LBorder+z,-yzero,0.2);
			glVertex3f(width-RBorder-z,-yzero,0.2);
			glColor4fv(col[i>0?0:1]);
			glVertex3f(width-RBorder,-yzero,0.2);
		}

		glEnd();

	}

	l = MIN(0,High);
	l = (int)(l/step)*step;
	if (l>=0)
		l= -step;

	for (double i=l;i>=Low-.001;i-=step) {
		float yzero=MapToY(i);

		glBegin(GL_LINE_STRIP);

		if (yzero>TBorder+2 && yzero<height-BBorder-2) {
			glColor4fv(col[0]);
			glVertex3f(LBorder,-yzero,0.2);
			glColor4fv(col[2]);
			glVertex3f(LBorder+z,-yzero,0.2);
			glVertex3f(width-RBorder-z,-yzero,0.2);
			glColor4fv(col[0]);
			glVertex3f(width-RBorder,-yzero,0.2);

		}

		glEnd();

	}

	//SetForeground(grey);



	glBegin(GL_LINE_STRIP);
	glColor4f(0,0,0,.5);
	glVertex3f(LBorder,-TBorder,0.25);
	glVertex3f(width-RBorder,-TBorder,0.25);
	glVertex3f(width-RBorder,-height+BBorder,0.25);
	glVertex3f(LBorder,-height+BBorder,0.25);
	glVertex3f(LBorder,-TBorder,0.25);
	glEnd();

	/*
	 glBegin(GL_LINE_STRIP);
	 glColor4f(0,0,0,.1);
	 glVertex3f(0,0,0.1);
	 glVertex3f(0,-height,0.1);
	 glVertex3f(width,-height,0.1);
	 glVertex3f(width,0,0.1);
	 glVertex3f(0,0,0.1);
	 glEnd();*/

}

float ToolPlot::MapFloatToX(double d) const {

	float ww = width-LBorder-RBorder;
	float w2 = LBorder;

	return (float)((d-currentView->xmin)/(currentView->xmax-currentView->xmin)*ww+.5f) + w2;
}

float ToolPlot::MapToY(double d) const {

	//int height = this->height - theight;

	return ((float)((High-d)*yscale)) + TBorder;
}

void ToolPlot::calcYscale()  {

	int height = this->height - theight;
	float hh = height-TBorder-BBorder;

	yscale = hh/(High-Low);
}


void ToolPlot::DrawString(int x, int y, const std::string &s) const {
	//glRasterPos3f(x,y,.1);
	Tool::font->printString(s.c_str(),x,y,.2);
}

double ToolPlot::getStep() const {

	double step = 30;
	if (currentView->xmax!=90 || currentView->xmin!=-90) {

		step = 25.*(currentView->xmax-currentView->xmin)/(width-RBorder-LBorder);
		step = roundStep(step);
	}
	else {
		if (width-RBorder-LBorder <100) step = 90;
		else if (width-RBorder-LBorder <150) step = 45;
		else if (width-RBorder-LBorder <200) step = 30;
		else if (width-RBorder-LBorder <400) step = 15;
		else if (width-RBorder-LBorder <700) step = 10;
		else/* if (width-RBorder-LBorder <800)*/ step = 5;


	}

	return step;
}

void ToolPlot::drawFrameText() const {

	int height = this->height - theight;

	//listbox->Draw(isActive,mouse);
	//DrawString(LBorder,-15,currentView->title);

	double step = getStep();



	for (double i=0;i<=currentView->xmax+.001;i+=step) {
		float xzero=MapFloatToX(i);


		if (xzero>=LBorder-2 && xzero<=width-RBorder+2) {
			std::string buf = ticktext(i);
			//   SetForeground(grey);
			//if (buf=="2")
			// int aaa = 0;

			int ww = Tool::font->StringLength(buf.c_str());

			//int pos = xzero-buf.length()*3+1;
			int pos = (int)xzero-ww/2;
			DrawString(pos,-height+BBorder-10-3,buf);
		}
	}

	for (double i=-step;i>=currentView->xmin-.001;i-=step) {
		float xzero=MapFloatToX(i);

		if (xzero>=LBorder-2 && xzero<=width-RBorder+2) {
			std::string buf = ticktext(i);
			// SetForeground(grey);
			//int pos = xzero-buf.length()*3+1;

			int ww = Tool::font->StringLength(buf.c_str());
			int pos = (int)xzero-ww/2;
			DrawString(pos,-height+BBorder-10-3, buf);
		}
	}


	step = 15*(High-Low)/(height-TBorder-BBorder);
	step = roundStep(step);

	double l = MAX(0,Low);
	l = (int)(l/step)*step;

	for (double i=l;i<=High+.001;i+=step) {
		float yzero=MapToY(i);

		if (yzero>=TBorder-2 && yzero<=height-BBorder+2) {
			//  SetForeground(grey);

			std::string buf = ticktext2(step,i);

			int ww = Tool::font->StringLength(buf.c_str());
			int pos = LBorder-ww-2;

			DrawString(pos,(int)-yzero-3,buf);
		}
	}

	l = MIN(0,High);
	l = (int)(l/step)*step;
	if (l>=0.)
		l= -step;

	for (double i=l;i>=Low-.001;i-=step) {
		float yzero=MapToY(i);

		if (yzero>=TBorder-2 && yzero<=height-BBorder+2) {

			std::string buf = ticktext2(step,i);

			int ww = Tool::font->StringLength(buf.c_str());
			int pos = LBorder-ww-2;
			//  SetForeground(grey);

			DrawString(pos,(int)-yzero-3,buf);
		}
	}
}


double roundHiLo(double diff) {
	double rd;
	if (diff<=.1)
		rd = .01;
	else if (diff<=.2)
		rd = .02;
	else if (diff<=.5)
		rd = .05;
	else if (diff<=1.)
		rd = .1;
	else if (diff<=2.)
		rd = .2;
	else if (diff<=5.)
		rd = .5;
	else if (diff<=10.)
		rd = 1.;
	else if (diff<=20.)
		rd = 2.;
	else if (diff<=50.)
		rd = 5.;
	else if (diff<=100.)
		rd = 10.;
	else if (diff<=200.)
		rd = 20.;
	else if (diff<=500.)
		rd = 50.;
	else if (diff<=1000.)
		rd = 100.;
	else if (diff<=2000.)
		rd = 200.;
	else /*if (diff<=5000.)*/
		rd = 500.;
	return rd;
}


void ToolPlot::newHighLow() {
	double rd, tmp;

	if (currentView) {
		currentView->GetMinMax(Low,High);
	}

	rd = roundHiLo(High-Low);
	tmp = ((int)(High/rd + .5))*rd;
	// printf("Hi: %f %f %f ",rd,0,tmp);
	if (High > tmp+.0001)
		High = tmp + rd;
	else
		High = tmp;
	tmp = ((int)(Low/rd - .5))*rd;
	// printf("Lo: %f %f %f\n",rd,0,tmp);
	if (Low >= tmp-.0001)
		Low = tmp;
	else
		Low = tmp - rd;

	if (High>10000)
		High=10000;
	if (Low<-10000)
		Low=-10000;
	if (ABS(High-Low)<.01) {
		High +=.1;
		Low -=.1;
	}

}


void ToolPlot::Apply() {

	std::vector<ToolItem*>::const_iterator it2 = items.begin();
	for (;it2!=items.end();it2++) {
		(*it2)->Apply();
	}

}

void ToolPlot::Discard() {

	std::vector<ToolItem*>::const_iterator it2 = items.begin();
	for (;it2!=items.end();it2++) {
		(*it2)->Discard();
	}

}

void ToolPlot::Draw(bool isActive, int HT, const Point<int> &mouse, bool border) {

	ToolWindow::Draw(isActive, HT, mouse, border);

	if (!currentView)
		return;

	// Projektion f�r inneres (mit Clip + Scroll)
	innerProj(HT);

	newHighLow();
	calcYscale();

	drawFrame();

	Point<int> m2 = mouse - Point<int>(x,y+theight-scrollcur);

	std::vector<ToolItem*>::const_iterator it2 = items.begin();
	for (;it2!=items.end();it2++) {
		(*it2)->Draw(/*activeItem == it2 &&*/ isActive, mouse);
	}

	drawFrameText();

	currentView->Draw(this);

}



bool ToolPlot::mousebutton(int mx, int my, int button, int type) {

	std::vector<ToolItem*>::iterator it2 = items.begin();
	for (;it2!=items.end();it2++) {
		if ((*it2)->TakesKey()) {
			if ((*it2)->contains(mx,my)) {
				//   hiliteItem = (*it2);

				(*it2)->MouseClick(button,type,Point<int>(mx,my));

				viewit = views.begin() + itemsel;

				currentView = *viewit;
				newHighLow();
				calcYscale();
				//activeItem = it2;
				return true;
			}
		}
	}


	if (button==1 && type == SDL_MOUSEBUTTONDOWN) {

		if (mx>LBorder && mx<width-RBorder && my>TBorder && my<height-BBorder) {
			viewit++;
			itemsel++;
			if (viewit==views.end()) {
				viewit = views.begin();
				itemsel = 0;
			}

			Discard();

			currentView = *viewit;

			newHighLow();
			calcYscale();
		}
	}
	/*
	 if (button==3 && type == SDL_MOUSEBUTTONDOWN) {

	  gl->iw->WriteTool(this);

	  cout << "KLICK" << endl;
	 }*/

	return true;
}


std::string ToolPlot::GetFilename() {

	if (currentView)
		return currentView->filename;

	return "test";

}

void ToolPlot::Reset() {

	clearpointer(items);
	clearpointer(views);
}



void ToolPlot::Finish() {

	std::vector<std::string> listentry;

	std::vector<ToolViewInterface*>::iterator it = views.begin();
	for (;it!=views.end();it++) {
		listentry.push_back((*it)->title);
	}

	ToolItem *listbox = new ToolList<int>(ToolRect(LBorder,2,140), itemsel, listentry);

	//listbox->Signal(plotclicklist);
	items.push_back(listbox);

	itemsel = 0;
	Discard();
	viewit = views.begin();
	currentView = *(viewit);
}

void ToolPlot::Resize() {

}
