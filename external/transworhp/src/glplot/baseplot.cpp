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

#include "baseplot.h"

#include "xopt_eps.h"


#include "../toolbase/tool.h"

#include "../core/Viewer.h"

#include "../base/vectortools.h"

#include <iostream>
#include <algorithm>

#include <GL/gl.h>



using namespace std;

namespace tw {

int BasePlot::LBorder=40;
int BasePlot::RBorder=7;
int BasePlot::TBorder=15;
int BasePlot::BBorder=7;
int BasePlot::minwidth=70;
int BasePlot::minheight=50;

int BasePlot::allplotnames=0;
int BasePlot::acsize=8;

//extern int onlyhigh;
BasePlot::BasePlot(int ind)
	: geom(100,100),
	  redflag(0),
	  extx(1),exty(1),
	  High(1.0), Low(0.0), High2(2.0), Low2(0.0),
	  index(ind),
	  mouseOver(0),
	  ctrldata(0),
	  comp(nullptr), comptime(nullptr),
	  dyncomp(nullptr), dyncomptime(nullptr),
	  dyncompn(nullptr),dyncompstep(nullptr),
	  scaledata(1.0),
	  zoomFactorX(0.0), zoomFactorY(0.0),
	  Mirror(false)
	   {

	calcYscale();

}


BasePlot::~BasePlot() {}
/*
void BasePlot::AddCompareCurve(double *cmp, int cmpstep) {

	comp = cmp;
	compstep = cmpstep;

}*/

void BasePlot::AddCompareCurve2(double *time, double *cmp, int cmpstep, int n) {

	comp = cmp;
	compstep = cmpstep;
	comptime = time;
	compn = n;

}

void BasePlot::AddDynCompareCurve(double *time, double *cmp, int *cmpstep, int *n) {

	dyncomp = cmp;
	dyncompstep = cmpstep;
	dyncomptime = time;
	dyncompn = n;

}

float jj[4][8][2] = {
	{
		{.5f,.5f}
	},
	{
		{.25f,.75f},{.75f,.25f}
	},

	{
		{.375f,.25f},{.125f,.75f},{.875f,.25f},{.625f,.75f}
	},

	{
		{.5625f,.4375f},{.0625f,.9375f},{.3125f,.6875f},{.6875f,.8125f},
		{.8125f,.1875f},{.9375f,.5625f},{.4375f,.0625f},{.1875f,.3125f}
	}
};


bool BasePlot::IsIcon() const {
	int width = GetWidth();
	int height = GetHeight();

	return (width==minwidth && height==minheight);
}

void BasePlot::Draw(DataStorage &ds, DataStorage &dstop, int ii) {

	Point<int> pos = GetPos();

	glPushMatrix();
	glTranslatef(pos.x,pos.y,0.f);

	if (IsIcon()) {
		drawIconFrame();
	} else {
		drawFrame();
	}

	//glEnable(GL_LINE_SMOOTH);
	//glDisable (GL_BLEND);
	glEnable (GL_BLEND);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

	/*glEnable (GL_LINE_SMOOTH);
	glEnable (GL_BLEND);
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	//glBlendFunc (GL_SRC_COLOR, GL_ONE_MINUS_SRC_COLOR);
	glHint (GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
	//glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	*/
//glLineWidth (1.5);


	/** Alternativ: Mit Multisampling -- aber leider nicht so glatt :-(
	glEnable(GL_MULTISAMPLE);
			glPushMatrix();
			glTranslatef(0,0,.1(i/(sz*2.)));
			drawData(ds, dstop);
			glPopMatrix();
	glDisable(GL_MULTISAMPLE);
	*/

	int sz = acsize;
	int jindex;

	switch (sz) {
	case 1:
		jindex=0;
		break;
	case 2:
		jindex=1;
		break;
	case 4:
		jindex=2;
		break;
	case 8:
	default:
		jindex=3;
		break;
	}

	for (int i=0; i<sz; i++) {
		glPushMatrix();
		glTranslatef(jj[jindex][i][0]-.5,jj[jindex][i][1]-.5,(i/(sz*2.)));
		drawData(ds, dstop,ii);
		glPopMatrix();
	}

	glDisable(GL_LINE_SMOOTH);
	/*glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	*/

	if (IsIcon()) {
		//drawIconFrame();
	} else {
		drawPlotName();
	}

	glPopMatrix();
}

void BasePlot::DrawMore(SDLFrame*, double*, double) {}

void BasePlot::drawText() const {

	if (!IsIcon()) {
		Point<int> pos = GetPos();
		glPushMatrix();
		glTranslatef(pos.x,pos.y,0.f);

		glColor4f(0.0f,0.0f,0.0f,1.0f);
		drawFrameText();

		glPopMatrix();
	}
}

void BasePlot::DrawText00(double, double) const {}
void BasePlot::drawData(DataStorage&, DataStorage&, int) {}
void BasePlot::SetMaxTime(DataStorage::TimeMode_e, double) {}
void BasePlot::SetMaxTime2(double) {}
void BasePlot::SetMinTime(DataStorage::TimeMode_e, double) {}

void BasePlot::SetUserControl(FunktionenzeigerC) {
		//ctrldata=1;
		//FuncC=c;
}

bool BasePlot::MouseInput(DataStorage&, int, const  Point<int>&) {
	return false;
}

bool BasePlot::MouseWheelInput(int scale, const  Point<int>& /*p*/) {
  
	if (zoomFactorX == 0.0 && zoomFactorY == 0.0) {
		//zoomFactorX = 25*(High2-Low2)/(GetWidth()-RBorder-LBorder);
		zoomFactorY = 15*(High-Low)/(GetHeight()-TBorder-BBorder);
	}

	// abhaengig von der Maus-Pos zoomen
	//double dev = static_cast<double>(GetHeight()-p.y)/GetHeight();

	// Y-Achse
	High += scale * zoomFactorY/2;
	Low -= scale * zoomFactorY/2;
	// X-Achse
	//High2 += scale * zoomFactorX/2;
	//Low2 -= scale * zoomFactorX/2;
	return true;
}

bool BasePlot::SpecialKeys(const Uint8*) {
	return false;
}

void BasePlot::TimerFunc(double) {}
void BasePlot::Info() {}

void BasePlot::SetEpsTitle(std::string buf) {
	epstitle = std::move(buf);
}
void BasePlot::SetScaleData(double sd) {
	scaledata = sd;
}
void BasePlot::ControlData(int i) {
	ctrldata = i;
}

void BasePlot::MoveTo(int currenttime, const Geometry &g, int time) {
	geom.SetNextGeometry(currenttime, g,time);
}
Point<int> BasePlot::GetPos() const {
	return geom.cur.pos;
}
int BasePlot::GetWidth() const {
	return geom.cur.width;
}
int BasePlot::GetHeight() const {
	return geom.cur.height;
}

int BasePlot::hasControlData() const {
	return ctrldata;
}

void BasePlot::matlab(std::ostream&) const {}

double BasePlot::roundStep(double step) const {

	double rd = 1.0;
	double base = 0.00001;

	for (int i=-2; i<20; i++) {

		if (step < base) {
			rd = base;
			break;
		} else if (step < 2*base) {
			rd = 2*base;
			break;
		} else if (step < 5*base) {
			rd = 5*base;
			break;
		}
		base *= 10;
	}

	return rd;
}

void BasePlot::calcYscale()  {

	const int height = GetHeight();

	float hh = height-TBorder-BBorder;

	if (IsIcon()) {
		hh = height;
	}

	yscale = hh/(High-Low);
}

float BasePlot::MapFloatToX(double d) const {
	const int width = GetWidth();

	float ww = static_cast<float>(width-LBorder-RBorder);
	float w2 = static_cast<float>(LBorder);
	if (IsIcon()) {
		ww = static_cast<float>(width);
		w2 = 0.0f;
	}
	return static_cast<float>( ((d-Low2)/(High2-Low2)*ww+.5f) + w2 );
}

double BasePlot::MapFromX(short x) const {
	const int width = GetWidth();

	float ww = static_cast<float>(width-LBorder-RBorder);
	float w2 = static_cast<float>(LBorder);
	if (IsIcon()) {
		ww = static_cast<float>(width);
		w2 = 0.0f;
	}
	return static_cast<double>(x-w2)/(ww+.5);
}

double BasePlot::MapFromY(short y) const {

	const int height = GetHeight();

	if (IsIcon()) {
		return High - (height - y)/yscale;
	} else {
		return High - (height - y - BBorder)/yscale;
	}
}



float BasePlot::MapTimeToX(DataStorage &ds, int i) const {

	const int width = GetWidth();

	float ww = static_cast<float>(width-LBorder-RBorder);
	float w2 = static_cast<float>(LBorder);
	if (IsIcon()) {
		ww = static_cast<float>(width);
		w2 = 0.0f;
	}

	//cout << ds.timemode << endl;

	if (ds.timemode==DataStorage::time_as_func) {
		float s = static_cast<float>(CallTimeFunc(i,ds.X.D,ds.UNKNOWN.D,ds.T.D,ds.X.ndis,ds.X.size,ds.UNKNOWN.size) / ds.getTF());
		//cout << "TIME " << i << " " << gl->ds.GetTF()<< " " << s << endl;
		return (s*ww+.5f)+w2;
	}

	else if (ds.timemode==DataStorage::time_is_time) {

		const float s = static_cast<float>( (ds.getData(Selector('t',0),i) - ds.MinT0)/(ds.MaxTF-ds.MinT0) );
		return (s*ww+.5f)+w2;
	}

	else if (ds.timemode==DataStorage::index_as_time) {
		//cout << "--- "  << i << " " << gl->ds.GetData(0,i)*gl->ds.GetTF() << endl;
		//return (float)(gl->ds.GetData(0,i)*ww+.5)+w2;
		int a = ds.getTF();
		//cout << "i " << i << " " << a << " " << i-1
		return static_cast<float>((i/(a-1.)*ww+.5)+w2);
	}

	else if (ds.timemode==DataStorage::dgl2_as_time) {

		double t = ds.getData(Selector('t',0),i)*ds.getData(Selector('x',DataStorage::timedgl),i);

		if (i >= DataStorage::nsplitdis) {

			int block = ((i-DataStorage::nsplitdis)/DataStorage::nsplitdis2);
			int b = DataStorage::nsplitindex[block]-1;
			double dt = ds.getData(Selector('t',0),b)*ds.getData(Selector('x',DataStorage::timedgl),b);

			/*			if ((i-DataStorage::nsplitdis)%DataStorage::nsplitdis2==0)
							absetzen=true;
						else
							absetzen=false;*/
			//std::cout << "XXX " << i << " "<< block << " " << b << " " << t << " " <<dt<< std::endl;
			t +=dt;

		}
		return static_cast<float>( (t/(ds.getTF())*ww+.5f)+w2 );

	}

	//cout << ds.start << " " << ds.tfgesamt << endl;
	//cout << "ffastime " << ds.floattime << " " <<  ds.tfgesamt<<" " << ds.timemode << endl;

	return static_cast<float>( ((ds.getData(Selector('t',0),i)*ds.getTF()/ds.MaxTF + ds.start)*ww+.5)+w2 );
}

float BasePlot::MapStopTimeToX(DataStorage &ds, DataStorage &dstop, int i) const {

	const int width = GetWidth();

	float ww = static_cast<float>(width-LBorder-RBorder);
	float w2 = static_cast<float>(LBorder);
	if (IsIcon()) {
		ww = static_cast<float>(width);
		w2 = 0.0f;
	}

	return ((dstop.getData(Selector('t',0),i)*dstop.getTF()+stopmark)/ds.MaxTF*ww+.5)+w2;
}


float BasePlot::MapToY(double d) const {
	//if (Mirror)
	//	return (short)((d-Low)*yscale)+TBorder;
	//else

	const int height = GetHeight();

	if (IsIcon()) {
		return static_cast<float>( height-(High-d)*yscale );
	} else {
		return static_cast<float>( height-(High-d)*yscale - BBorder );
	}
}

void BasePlot::drawFrame() const {

	const int width = GetWidth();
	const int height = GetHeight();

	// Hintergrund
	const std::vector<GLfloat> verticesBackground = {
		static_cast<GLfloat>(LBorder), static_cast<GLfloat>(TBorder), 0.0f,
		static_cast<GLfloat>(LBorder), static_cast<GLfloat>(height-BBorder), 0.0f,
		static_cast<GLfloat>(width-RBorder), static_cast<GLfloat>(TBorder), 0.0f,
		static_cast<GLfloat>(width-RBorder), static_cast<GLfloat>(height-BBorder), 0.0f
	};

	// Rahmen/Begrenzung direkt am Gitter
	const std::vector<GLfloat> verticesInnerBorder = {
		static_cast<GLfloat>(LBorder),		static_cast<GLfloat>(TBorder),		0.1f,
		static_cast<GLfloat>(width-RBorder),	static_cast<GLfloat>(TBorder),		0.1f,
		static_cast<GLfloat>(width-RBorder),	static_cast<GLfloat>(height-BBorder),	0.1f,
		static_cast<GLfloat>(LBorder),		static_cast<GLfloat>(height-BBorder),	0.1f,
	};

	// Rahmen um das Plot Fenster
	const std::vector<GLfloat> verticesFrameBorder = {
		0.0f,				0.0f,				0.1f,
		0.0f,				static_cast<GLfloat>(height),	0.1f,
		static_cast<GLfloat>(width),	static_cast<GLfloat>(height),	0.1f,
		static_cast<GLfloat>(width),	0.0f,				0.1f
	};

	// Gitterlinien
	std::vector<GLfloat> verticesGrid;
	verticesGrid.reserve(256); //heuristisch..

	double step = roundStep( 25*(High2-Low2)/(width-RBorder-LBorder) );

	// vertikale Linien
	for (double i = Low2-0.001; i <= High2+.001; i+=step) {

		const float xzero = MapFloatToX(i);

		//damit nicht aus dem (Plot-)Fenster gemalt wird
		if (xzero > LBorder+2 && xzero < width-RBorder-2) {
			verticesGrid.push_back(xzero);
			verticesGrid.push_back(TBorder);
			verticesGrid.push_back(0.1f);

			verticesGrid.push_back(xzero);
			verticesGrid.push_back(height-BBorder);
			verticesGrid.push_back(0.1f);
		}
	}

	step = roundStep( 15*(High-Low)/(height-TBorder-BBorder) );

	// horizontale Linien
	for (double i = Low-.001; i <= High+.001; i+=step) {

		const float yzero = MapToY(i);

		//damit nicht aus dem (Plot-)Fenster gemalt wird
		if (yzero > TBorder+2 && yzero < height-BBorder-2) {
			verticesGrid.push_back(LBorder);
			verticesGrid.push_back(yzero);
			verticesGrid.push_back(0.1f);

			verticesGrid.push_back(width-RBorder);
			verticesGrid.push_back(yzero);
			verticesGrid.push_back(0.1f);
		}
	}

	const float alpha = GetHue()*.2f+.3f;

	glBindBuffer(GL_ARRAY_BUFFER, 0);

	glEnableClientState(GL_VERTEX_ARRAY);

	glColor4f(1.0f, 1.0f, 1.0f, alpha);
	glVertexPointer(3, GL_FLOAT, 0, verticesBackground.data());
	glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);

	glColor4f(0.0f, 0.0f, 0.0f, 0.2f);
	glVertexPointer(3, GL_FLOAT, 0, verticesGrid.data());
	glDrawArrays(GL_LINES, 0, verticesGrid.size()/3);

	glColor4f(0.0f, 0.0f, 0.0f, 0.5f);
	glVertexPointer(3, GL_FLOAT, 0, verticesInnerBorder.data());
	glDrawArrays(GL_LINE_LOOP, 0, 4);

	glColor4f(0.0f, 0.0f, 0.0f, 0.1f);
	glVertexPointer(3, GL_FLOAT, 0, verticesFrameBorder.data());
	glDrawArrays(GL_LINE_LOOP, 0, 4);

	glDisableClientState(GL_VERTEX_ARRAY);
}

float BasePlot::GetHue() const {
	float th = (mouseOver-20)/30.f;
	if (th > 1.0f) {
		th = 1.0f;
	} else if (th < 0.0f) {
		th = 0.0f;
	}
	return th;
}

void BasePlot::drawPlotName() const {

	const int width = GetWidth();
	const int height = GetHeight();

	if (mouseOver > 20 || allplotnames) {

		float th = GetHue();

		if (allplotnames) {
			th = 1.0f;
		}

		if (epstitle.size()) {
			const int sl = Tool::font->StringLength(epstitle.c_str());

			const float ll = static_cast<float>( (width-LBorder-RBorder-sl)/2+LBorder );
			const float rr = static_cast<float>( ll+sl );
			const float tt = static_cast<float>( height-BBorder-10 );
			const float bb = static_cast<float>( height-BBorder+0 );

			const float h = .9f;

			const std::vector<GLfloat> verticesBackground = {
				ll-3.0f,	tt-3.0f,	h,
				rr+3.0f,	tt-3.0f,	h,
				ll-3.0f,	bb+3.0f,	h,
				rr+3.0f,	bb+3.0f,	h

			};

			const std::vector<GLfloat> colorBackground = {
				1.0f,	236.0f/255.0f,	139.0f/255.0f,	.6f*th,
				1.0f,	236.0f/255.0f,	139.0f/255.0f,	.6f*th,
				1.0f,	236.0f/255.0f,	139.0f/255.0f,	.9f*th,
				1.0f,	236.0f/255.0f,	139.0f/255.0f,	.9f*th
			};

			const std::vector<GLfloat> verticesBorder = {
				ll-3.0f,	tt-3.0f,	h+.1f,
				rr+3.0f,	tt-3.0f,	h+.1f,
				rr+3.0f,	bb+3.0f,	h+.1f,
				ll-3.0f,	bb+3.0f,	h+.1f,
			};

			const std::vector<GLfloat> colorBorder = {
				0.0f,	0.0f,	0.0f,	0.7f*th,
				0.0f,	0.0f,	0.0f,	0.7f*th,
				0.0f,	0.0f,	0.0f,	0.7f*th,
				0.0f,	0.0f,	0.0f,	0.7f*th
			};


			glBindBuffer(GL_ARRAY_BUFFER, 0);

			glEnableClientState(GL_VERTEX_ARRAY);
			glEnableClientState(GL_COLOR_ARRAY);

			glColorPointer(4, GL_FLOAT, 0, colorBackground.data());
			glVertexPointer(3, GL_FLOAT, 0, verticesBackground.data());
			glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);

			glColorPointer(4, GL_FLOAT, 0, colorBorder.data());
			glVertexPointer(3, GL_FLOAT, 0, verticesBorder.data());
			glDrawArrays(GL_LINE_LOOP, 0, 4);

			glDisableClientState(GL_COLOR_ARRAY);
			glDisableClientState(GL_VERTEX_ARRAY);

			glColor4f(139.0f/255.0f,129.0f/255.0f,76.0f/255.0f,th);
			Tool::font->printString(epstitle.c_str(),ll+1.0f,tt,h+.1f);
		}
	}
}


void BasePlot::drawIconFrame() const {

	const int width = GetWidth();
	const int height = GetHeight();

	// Hintergrund
	const std::vector<GLfloat> verticesBackground = {
		0.0f, 				0.0f,				0.0f,
		0.0f, 				static_cast<GLfloat>(height), 	0.0f,
		static_cast<GLfloat>(width), 	0.0f, 				0.0f,
		static_cast<GLfloat>(width), 	static_cast<GLfloat>(height),	0.0f
	};

	// Rahmen um das Plot Fenster
	const std::vector<GLfloat> verticesFrameBorder = {
		0.0f,				0.0f,				0.1f,
		0.0f,				static_cast<GLfloat>(height),	0.1f,
		static_cast<GLfloat>(width),	static_cast<GLfloat>(height),	0.1f,
		static_cast<GLfloat>(width),	0.0f,				0.1f
	};

	// Gitterlinien
	std::vector<GLfloat> verticesGrid;
	verticesGrid.reserve(64); //heuristisch..

	double step = roundStep( 30*(High2-Low2)/width );

	// vertikale Linien
	for (double i = Low2-0.001; i <= High2+.001; i+=step) {

		const float xzero = MapFloatToX(i);

		//damit nicht aus dem (Plot-)Fenster gemalt wird
		if (xzero > 1 && xzero < width-1) {
			verticesGrid.push_back(xzero);
			verticesGrid.push_back(0.0f);
			verticesGrid.push_back(0.1f);

			verticesGrid.push_back(xzero);
			verticesGrid.push_back(height);
			verticesGrid.push_back(0.1f);
		}
	}

	step = roundStep( 20*(High-Low)/height );

	// horizontale Linien
	for (double i = Low-.001; i <= High+.001; i+=step) {

		const float yzero = MapToY(i);

		//damit nicht aus dem (Plot-)Fenster gemalt wird
		if (yzero > 1 && yzero < height-1) {
			verticesGrid.push_back(0.0f);
			verticesGrid.push_back(yzero);
			verticesGrid.push_back(0.1f);

			verticesGrid.push_back(width);
			verticesGrid.push_back(yzero);
			verticesGrid.push_back(0.1f);
		}
	}

	const float alpha = GetHue()*.2f+.3f;

	glBindBuffer(GL_ARRAY_BUFFER, 0);

	glEnableClientState(GL_VERTEX_ARRAY);

	glColor4f(1.0f, 1.0f, 1.0f, alpha);
	glVertexPointer(3, GL_FLOAT, 0, verticesBackground.data());
	glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);

	glColor4f(0.0f, 0.0f, 0.0f, 0.2f);
	glVertexPointer(3, GL_FLOAT, 0, verticesGrid.data());
	glDrawArrays(GL_LINES, 0, verticesGrid.size()/3);

	glColor4f(0.0f, 0.0f, 0.0f, 0.1f);
	glVertexPointer(3, GL_FLOAT, 0, verticesFrameBorder.data());
	glDrawArrays(GL_LINE_LOOP, 0, 4);

	glDisableClientState(GL_VERTEX_ARRAY);
}



/** Check if a double (nearly) ends in .000 */
inline bool is_round(double value) {
	double v = value-(int)value;
	if (v<0) {
		v=-v;
	}
	return (v<1e-8);
}


/** Create Tick text (x) */
std::string ticktext(double step, double i) {
	std::stringstream s;

	if (step>1000) {

		if (i==0) {
			s << int( i );
		} else {
			//s.setf(ios::left, ios::adjustfield);
			s.setf(std::ios::scientific, std::ios::floatfield);

			s << std::setprecision(1) << std::setw(2) ;
			s << ( i );
		}
	} else {


		if (is_round(i)) {
			s << int( i );
		} else if (is_round(i*10)) {
			s << std::setw(1) << i;
		} else {
			s << std::setw(2) << i;
		}
	}
	return s.str();
}

/** Create Tick text (y) */
std::string ticktext2(double step, double i) {

	std::stringstream s;

	if (step >= 10000.0) {
		//s.setf(ios::left, ios::adjustfield);
		s.setf(std::ios::scientific, std::ios::floatfield);
		const int ii = static_cast<int>(i);
		s << std::setprecision(2) << std::setw(10) << static_cast<double>(ii);
	} else if (step >= 1.0) {
		s << std::setw(5) << static_cast<int>(i);
	} else if (step < 0.0001) {
		s << std::setiosflags(std::ios::showpoint | std::ios::right | std::ios::fixed)
		  << std::setprecision(5) << std::setw(7) << i;
	} else if (step < 0.001) {
		s << std::setiosflags(std::ios::showpoint | std::ios::right | std::ios::fixed)
		  << std::setprecision(4) << std::setw(6) << i;
	} else if (step < 0.01) {
		s << std::setiosflags(std::ios::showpoint | std::ios::right | std::ios::fixed)
		  << std::setprecision(3) << std::setw(5) << i;
	} else if (step < 0.1) {
		s << std::setiosflags(std::ios::showpoint | std::ios::right | std::ios::fixed)
		  << std::setprecision(2) << std::setw(5) << i;
	} else {
		s << setiosflags(std::ios::showpoint | std::ios::right | std::ios::fixed)
		  << std::setprecision(1) << std::setw(5) << i;
	}
	return s.str();
}


void BasePlot::drawFrameText() const {

	const int width = GetWidth();
	const int height = GetHeight();


	double step = 25*(High2-Low2)/(width-RBorder-LBorder);

	step = roundStep(step);


	for (double i = 0; i <= High2+.001; i+=step) {
		const float xzero = MapFloatToX(i);


		if (xzero >= LBorder-2 && xzero <= width-RBorder+2) {
			const std::string buf = ticktext(step,i);
			//			SetForeground(grey);
			//			if (buf=="2")
			//				int aaa = 0;

			const int ww = Tool::font->StringLength(buf.c_str());

			//int pos = xzero-buf.length()*3+1;
			const int pos = xzero-ww/2;
			DrawString(pos,TBorder-10-3,buf);
		}
	}

	for (double i = -step; i >= Low2-.001; i-=step) {
		const float xzero = MapFloatToX(i);

		if (xzero >= LBorder-2 && xzero <= width-RBorder+2) {
			const std::string buf = ticktext(step,i);
			//	SetForeground(grey);
			//int pos = xzero-buf.length()*3+1;

			const int ww = Tool::font->StringLength(buf.c_str());
			const int pos = xzero-ww/2;
			DrawString(pos,TBorder-10-3, buf);
		}
	}


	step = 15*(High-Low)/(height-TBorder-BBorder);
	step = roundStep(step);

	double k = std::max(0., Low);
	k = (int)(k/step)*step;

	for (double i = k; i <= High+.001; i+=step) {
		const float yzero = MapToY(i);

		if (yzero >= TBorder-2 && yzero <= height-BBorder+2) {
			//		SetForeground(grey);

			const std::string buf = ticktext2(step,i);

			const int ww = Tool::font->StringLength(buf.c_str());
			const int pos = LBorder-ww-2;

			DrawString(pos,yzero-5,buf);
		}
	}

	k = min(0., High);
	k = (int)(k/step)*step;
	if (k>=0.)
		k= -step;


	for (double i = k; i >= Low-.001; i-=step) {
		const float yzero=MapToY(i);

		if (yzero >= TBorder-2 && yzero <= height-BBorder+2) {

			const std::string buf = ticktext2(step,i);

			const int ww = Tool::font->StringLength(buf.c_str());
			const int pos = LBorder-ww-2;
			//		SetForeground(grey);

			DrawString(pos,yzero-5,buf);
		}
	}

}

void BasePlot::DrawString(int x, int y, const std::string &s) const {
	//glRasterPos3f(x,y,.1);
	Tool::font->printString(s.c_str(),x,y,.1);
}


bool BasePlot::contains(const Point<int> &p)  {

	return p.isin(GetPos(),GetPos()+Point<int>(GetWidth(),GetHeight()));
}

/*
1: Mitte
2: L
4: R
8: T
16: B
*/
int BasePlot::containsArea(const Point<int> &p)  {

	int ret = 0;
	if (contains(p)) {

		const int width = GetWidth();
		const int height = GetHeight();
		const Point<int> pos = GetPos();

		if (IsIcon()) {
			if (p.isin(pos, pos+Point<int>(2,height)))
				ret +=2;
			if (p.isin(pos+Point<int>(width-2,0), pos+Point<int>(width,height)))
				ret +=4;
			if (p.isin(pos+Point<int>(0,height-2), pos+Point<int>(width,height)))
				ret +=8;
			if (p.isin(pos, pos+Point<int>(width,2)))
				ret +=16;

			if (ret==0 && contains(p))
				ret+=1;
		} else {

			if (p.isin(pos, pos+Point<int>(LBorder,height)))
				ret +=2;
			if (p.isin(pos+Point<int>(GetWidth()-RBorder,0), pos+Point<int>(width,height)))
				ret +=4;
			if (p.isin(pos+Point<int>(0,height-BBorder), pos+Point<int>(width,height)))
				ret +=8;
			if (p.isin(pos, pos+Point<int>(width,TBorder)))
				ret +=16;

			if (ret==0 && contains(p))
				ret+=1;
		}
	}

	return ret;
}



void BasePlot::DrawDynCompareCurve(DataStorage &ds) const {

	if (!dyncomp) {
		return;
	}

	if (*dyncompstep == 0) {

		const int cury = MapToY(dyncomp[0]);

		SetColor(TWcolor::Green);
		glBegin(GL_LINE_STRIP);
		glVertex3f(0.0f, cury,0.7f);
		glVertex3f(LBorder, cury,0.7f);
		glEnd();
		glBegin(GL_LINE_STRIP);
		glVertex3f(GetWidth()-RBorder, cury,0.7f);
		glVertex3f(GetWidth(), cury,0.7f);

		glEnd();

	} else if (dyncomptime == nullptr) {
		SetColor(TWcolor::Green);

		glBegin(GL_LINE_STRIP);
		for (int i = 0; i < *dyncompn; i++) {

			const int cury = MapToY(dyncomp[*dyncompstep*i]);
			const int curx = MapTimeToX(ds,i);

			glVertex3f(curx, cury,0.7f);
		}
		glEnd();

	} else {
		SetColor(TWcolor::Green);
		//std::cout << "COMPAREC" << comptime[40] <<  std::endl;

		glBegin(GL_LINE_STRIP);
		for (int i = 0; i < *dyncompn; i++) {

			//std::cout << "COMPAREC" << dyncomp[*dyncompstep*i]
			//<< " " <<  dyncomptime[i] <<  std::endl;
			const int cury = MapToY(dyncomp[*dyncompstep*i]);
			const int curx = MapFloatToX(dyncomptime[i]);

			glVertex3f(curx, cury,0.7f);
		}
		glEnd();
	}
}


void BasePlot::DrawCompareCurve(DataStorage &ds) const {

	if (!comp) {
		return;
	}

	// MR:was macht dieser Fall? Linie links und rechts vom Plot?
	if (compstep == 0) {

		const int cury = MapToY(comp[0]);

		SetColor(TWcolor::Grey);
		glBegin(GL_LINE_STRIP);
		glVertex3f(0.0f, static_cast<float>(cury), 0.5f);
		glVertex3f(static_cast<float>(LBorder), static_cast<float>(cury), 0.5f);
		glEnd();

		glBegin(GL_LINE_STRIP);
		glVertex3f(static_cast<float>(GetWidth()-RBorder), static_cast<float>(cury), 0.5f);
		glVertex3f(static_cast<float>(GetWidth()), static_cast<float>(cury), 0.5f);
		glEnd();

	} else if (comptime == nullptr) {
		SetColor(TWcolor::Grey);

		glBegin(GL_LINE_STRIP);
		for (int i = 0; i < compn; i++) {

			const int cury = MapToY(comp[compstep*i]);
			const int curx = MapTimeToX(ds,i);

			glVertex3f(static_cast<float>(curx), static_cast<float>(cury), 0.5f);
		}
		glEnd();

	} else {
		SetColor(TWcolor::Grey);
		//std::cout << "COMPAREC" << comptime[40] <<  std::endl;

		glBegin(GL_LINE_STRIP);
		for (int i = 0; i < compn; i++) {

			//	std::cout << "COMPAREC" << comp[compstep*i]
			//	<< " " <<  comptime[i] <<  std::endl;
			const int cury = MapToY(comp[compstep*i]);
			const int curx = MapFloatToX(comptime[i]);

			glVertex3f(static_cast<float>(curx), static_cast<float>(cury), 0.5f);
		}
		glEnd();
	}
}


/** Round to 1,2 or 5 */
double roundHiLo(double diff) {
	double rd = 0.0;

	double base = 0.00001;

	for (int i = -2; i < 20; i++) {

		if (diff <= base) {
			rd = base * .1;
			break;
		} else if (diff <= 2.*base) {
			rd = base * .2;
			break;
		} else if (diff <= 5.*base) {
			rd = base * .5;
			break;
		}
		base *= 10;
	}

	return rd;
}

void BasePlot::newHighLow(DataStorage &ds) {

	int ret = rawHighLow(ds);
	if (ret == -1) return;

	if (comp && compstep) {

		for (int i = 0; i < compn; i++) {
			double cury = comp[compstep*i];
			if (High < cury)
				High = cury;
			if (Low > cury)
				Low = cury;
		}
	}

	if (dyncomp && *dyncompstep) {

		for (int i = 0; i < *dyncompn; i++) {
			double cury = dyncomp[*dyncompstep*i];
			if (High < cury)
				High = cury;
			if (Low > cury)
				Low = cury;
		}
	}


	double rd = roundHiLo((High-Low)/2.1);


	//std::cout << setw(5) << index << ":  " << setw(15) << Low << setw(15) << High <<  setw(15) << rd  << endl;

	double aa = (High/rd + .5);

	double tmp;
	if (fabs(aa) > 1e10) {
		tmp = (int)(High*1.1);
	} else {
		tmp = ((int)aa)*rd ;
	}
	//	printf("Hi: %f %f %f ",rd,0,tmp);

	//if (tmp<0) {
	//	tmp*=-1;
	//	cout << "Bad tmp" << High<< " " << rd << " " <<High/rd+.5<<  " " << ((long)(High/rd + .5))*rd << endl;
	//}

	if (High > tmp+.0001) {
		High = tmp + rd;
	} else {
		High = tmp;
	}

	aa = (Low/rd - .5);
	if (fabs(aa)>1e10) {
		tmp = (int)(Low*0.9);
	} else {
		tmp = ((int)aa)*rd ;
	}

	//tmp = ((int)(Low/rd - .5))*rd;
	//	printf("Lo: %f %f %f\n",rd,0,tmp);

	//std::cout << index << ":  " << Low << " " << tmp  << " " << Low -tmp<< std::endl;

	if (Low >= tmp-.0001)
		Low = tmp;
	else
		Low = tmp - rd;


	//std::cout << setw(10) << Low << setw(10) << High << std::endl;

	if (High>1e20)
		High=1e20;
	if (Low<-1e20)
		Low=-1e20;
	if (fabs(High-Low)<.0001) {
		High +=.001;
		Low -=.001;
	}

	//std::cout << index << ":  " << Low << " " << High << " " << rd << std::endl;

	//printf("%f %f\n",High,Low);
	if (!ret)
		return;


	rd = roundHiLo(High2-Low2);
	tmp = ((int)(High2/rd + .5))*rd;
	//printf("Hi: %f %f %f ",rd,High2,tmp);
	if (High2 > tmp+.0001)
		High2 = tmp + rd;
	else
		High2 = tmp;
	tmp = ((int)(Low2/rd - .5))*rd;
	//printf("Lo: %f %f %f\n",rd,Low2,tmp);
	if (Low2 >= tmp-.0001)
		Low2 = tmp;
	else
		Low2 = tmp - rd;

	//std::cout << index << "::  " << Low2 << " " << High2 << std::endl;

}



void BasePlot::SetColor(TWcolor c) const {

	//glColor4f(0.f,0.f,0.f,0.f);
	//return;

	//std::cout << "AA " <<  gl->acsize << std::endl;

	const int sz = acsize*.8;

	switch (c) {
	case TWcolor::Black:
		glColor4f(0.0f,0.0f,0.0f,1.0f/sz);
		break;
	case TWcolor::White:
		glColor4f(0.0f,0.0f,0.8f,1.0f/sz);
		break;
	case TWcolor::Red:
		glColor4f(0.8f,0.0f,0.0f,1.0f/sz);
		break;
	case TWcolor::Rose:
		glColor4f(0.8f,0.4f,0.4f,1.0f/sz);
		break;
	case TWcolor::Grey:
		glColor4f(0.5f,0.5f,0.5f,1.0f/sz);
		break;
	case TWcolor::Green:
		glColor4f(0.0f,0.8f,0.0f,1.0f/sz);
		break;
	case TWcolor::Cyan:
		glColor4f(0.4f,0.4f,1.0f,1.0f/sz);
		break;
	case TWcolor::Hue: {
		const float hue = .5f; // gl->hue;

		int sz = acsize*.8;
		float a[3] = {0.,100./255.,0.};
		float b[3] = {0.,.8,0.};
		float c[4];
		for (int i=0; i<3; i++) {
			c[i] = a[i] * (1.0f-hue) + b[i] * (hue);
		}
		c[3] = 1.0f/sz;
		glColor4fv(c);
		}
		break;
	case TWcolor::Hue2: {
		const float hue = .5f; // gl->hue;
		int sz = acsize*.8;
		float a[3] = {100./255.,0.,0.};
		float b[3] = {.8,0.,0.};
		float c[4];
		for (int i=0; i<3; i++) {
			c[i] = a[i] * (1.0f-hue) + b[i] * (hue);
		}
		c[3] = 1.0f/sz;
		glColor4fv(c);
		}
		break;
	}
}


void BasePlot::Timer(Viewer *v) {

	if (geom.Update(v->change, v->currentTime)) {
		calcYscale();
	}

	if (v->mouseOver == this) {
		mouseOver+=5;
		if (mouseOver>50)
			mouseOver=50;
	} else {
		if (mouseOver>0)
			mouseOver--;
	}
}

BasePlot::Geometry::Geometry(int w, int h)
	: width(w), height(h) {}

BasePlot::DynGeometry::DynGeometry(int w, int h)
	: last(w,h), next(w,h), cur(w,h) {
	lasttime = 0;
	nexttime = 0;

	//cout << "BasePlot:DynG" << nexttime <<  endl;
}


bool BasePlot::DynGeometry::Update(int &change, int currenttime) {

	if (nexttime == 0) {
		return false;
	}

	if (currenttime < nexttime + 100) {
		change = 100;
	}

	//std::cout << gl->Time<< endl;
	if (currenttime > nexttime) {
		cur = next;
		return true;
	} else {
		const double a = static_cast<double>(currenttime-lasttime)/(nexttime-lasttime);
		cur.width = last.width+a*(next.width-last.width);
		cur.height = last.height+a*(next.height-last.height);
		cur.pos.x = last.pos.x+a*(next.pos.x-last.pos.x);
		cur.pos.y = last.pos.y+a*(next.pos.y-last.pos.y);
		return true;
	}
}

void BasePlot::DynGeometry::SetNextGeometry(int currenttime, const BasePlot::Geometry &g, int time) {
	last = cur;
	next = g;
	lasttime = currenttime;
	nexttime = lasttime+time;

	//std::cout << "SetN"<< lasttime << " " << nexttime << " " << g.width<< endl;
}

}
