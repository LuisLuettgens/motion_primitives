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
#include "toolview.h"
#include "toolplot.h"

#include <iomanip>
#include "GL/gl.h"
#include "conversion.h"
#include "../base/vectortools.h"

using namespace std;


#define MAX(A,B) ((A)>(B)?(A):(B))
#define MIN(A,B) ((A)<(B)?(A):(B))
#define ABS(A) ((A)>=0?(A):-(A))


ToolViewInterface::ToolViewInterface(const std::string &t, const std::string &f, int l2,int h2,
									 int nx, int ny)
		: xmin(l2), xmax(h2), ndis(nx), ndgl(ny), title(t), filename(f) {}


void ToolViewInterface::GetMinMax(double &Low, double &High) {
	for (int i=0;i<ndis;i++) {
		double vy = Get(i,ndgl-1);

		if (i==0 || vy < Low)
			Low = vy;
		if (i==0 || vy > High)
			High = vy;

	}
}

/*  --------------------------------------------- */

template <class T>
ToolViewLines<T>::ToolViewLines(const std::string &t, const std::string &f, int l2,int h2, T *vy,
								int nx)
	: ToolViewInterface(t, f,l2,h2,nx,1), valy(vy) {}

template <class T>
double ToolViewLines<T>::Get(int x, int y) {
	return valy[x*ndgl+y];
}


template <class T>
void ToolViewLines<T>::Draw(ToolPlot *p) {

	glBegin(GL_LINES);

	for (int i=0;i<ndis;i++) {
		float x = p->MapFloatToX(i*(xmax-xmin)/(ndis-1.)+xmin);
		float y = p->MapToY(Get(i,0));
		float y0 = p->MapToY(0);

		glColor4f(0,0,.8,.5);
		glVertex3f(x,-y,0.3);
		glVertex3f(x,-y0,0.3);
	}

	glEnd();
}

template <class T>
void ToolViewLines<T>::Text(std::ofstream &of) {

	for (int i=0;i<ndis;i++) {
		of << setw(20) << i*(xmax-xmin)
		/(ndis-1.)+xmin;
		of << setw(20) << Get(i,0);
		of << endl;
	}
}



/*  --------------------------------------------- */




template <class T>
ToolViewBars<T>::ToolViewBars(const std::string &t, const std::string &f, int l2,int h2, T *vy,
							  int nx, int nor)
		: ToolViewInterface(t, f,l2,h2,nx,1), valy(vy), norm(nor) {}


template <class T>
double ToolViewBars<T>::Get(int x, int y) {
	return valy[x*ndgl+y];
}

template <class T>
void ToolViewBars<T>::GetMinMax(double &Low, double &High) {
	for (int i=0;i<ndis;i++) {
		double vy = Get(i,ndgl-1);

		if (norm)
			vy *= 100./norm;

		if (vy < Low || i==0)
			Low = vy;
		if (vy > High || i==0)
			High = vy;
	}
}


template <class T>
void ToolViewBars<T>::Draw(ToolPlot *p) {

	glBegin(GL_QUADS);

	for (int i=0;i<ndis;i++) {
		float x0 = p->MapFloatToX(i-.5);
		float x1 = p->MapFloatToX(i+.5)-1;

		double nn = Get(i,0);
		if (norm)
			nn *= 100./norm;

		float y = p->MapToY(nn);
		float y0 = p->MapToY(0);

		glColor4f(0,0,.8,.5);
		glVertex3f(x0,-y,0.3);
		glVertex3f(x0,-y0,0.3);
		glVertex3f(x1,-y0,0.3);
		glVertex3f(x1,-y,0.3);

	}

	glEnd();
}

template <class T>
void ToolViewBars<T>::Text(std::ofstream &of) {

	for (int i=0;i<ndis;i++) {

		double nn = Get(i,0);
		if (norm)
			nn *= 100./norm;

		of << setw(20) << i;
		of << setw(20) << nn;
		of << endl;
	}

}


/*  --------------------------------------------- */



ToolViewGraph::ToolViewGraph(const std::string &t, const std::string &f, int l2,int h2,
							 std::vector<Point<double> > *vy, const std::stringstream &str)
	: ToolViewInterface(t, f,l2,h2,0,0), valy(vy), text(str.str()) {}



void ToolViewGraph::Draw(ToolPlot *p) {

	if (this->valy->size()<2) return;


	std::vector<Point<double> >::iterator it = this->valy->begin();

	glColor4f(.8,0,.0,.5);

	glBegin(GL_LINE_STRIP);
	int x = (int)p->MapFloatToX(it->x);
	int y = (int)p->MapToY(it->y);
	glVertex3f(x,-y,0.3);
	it++;
	x = (int)p->MapFloatToX(it->x);
	y = (int)p->MapToY(it->y);
	glVertex3f(x,-y,0.3);
	glEnd();
	
	char buf[20];
	sprintf(buf,"%.3f",it->y);
		
	glColor4f(.8,0,.0,.9);
	int a = Tool::font->StringLength(buf);
	Tool::font->printString(buf,x-a,-y+1,.35);
	
	
	it++;

	
	glColor4f(0,0,.8,.5);

	glEnable(GL_LINE_SMOOTH);
	glLineWidth(2);

	glBegin(GL_LINE_STRIP);

	for (;it!=this->valy->end();it++) {

		x = (int)p->MapFloatToX(it->x);
		y = (int)p->MapToY(it->y);
		glVertex3f(x,-y,0.4);
	}

	glEnd();
	glLineWidth(1);

	glDisable(GL_LINE_SMOOTH);

}


void ToolViewGraph::Text(std::ofstream &of) {

	of << text;
}

void ToolViewGraph::GetMinMax(double &Low, double &High) {

	if (this->valy->size()<2) return;

	std::vector<Point<double> >::iterator it = this->valy->begin();

	if (it!=valy->end()) {
		High = it->y;
		Low = it->y;
		it++;
	}

	for (;it!=this->valy->end();it++) {

		if (High<it->y) High = it->y;
		if (Low>it->y) Low = it->y;

	}

}


/*  --------------------------------------------- */



template <class T, int N>
ToolViewMinMax<T,N>::ToolViewMinMax(const std::string &t, const std::string &f,
									int l2,int h2,MinMax<T,N>  *vy)
	: ToolViewInterface(t, f,l2,h2,N,3), valy(vy) {}


template <class T, int N>
double ToolViewMinMax<T,N>::Get(int x, int y) {
	return valy->Get(x,y);
}


template <class T, int N>
void ToolViewMinMax<T,N>::Draw(ToolPlot *p) {

	glBegin(GL_LINES);

	for (int i=0;i<ndis;i++) {
		float x = p->MapFloatToX(i*(xmax-xmin)/(ndis-1.)+xmin);
		float y = p->MapToY(Get(i,1));
		float y0 = p->MapToY(Get(i,2));

		glColor4f(0,0,.8,.5);
		glVertex3f(x,-y,0.3);
		glVertex3f(x,-y0,0.3);
	}

	glEnd();

	glPointSize(3);
	glEnable(GL_POINT_SMOOTH);
	glBegin(GL_POINTS);


	for (int i=0;i<ndis;i++) {
		float x = p->MapFloatToX(i*(xmax-xmin)/(ndis-1.)+xmin);
		float y = p->MapToY(Get(i,0));
		glColor4f(0,0,.2,.3);
		glVertex3f(x,-y,0.4);

	}

	glEnd();

	glPointSize(1);
	glDisable(GL_POINT_SMOOTH);
}

template <class T, int N>
void ToolViewMinMax<T,N>::Text(std::ofstream &of) {

	for (int i=0;i<ndis;i++) {

		of << setw(20) << i*(xmax-xmin)
		/(ndis-1.)+xmin;
		of << setw(20) << Get(i,1);
		of << setw(20) << Get(i,2);
		of << setw(20) << Get(i,0);
		of << endl;
	}

}


/*  --------------------------------------------- */




template <class T>
ToolViewColor<T>::ToolViewColor(const std::string &t, const std::string &f, int l2,int h2, T *vy,
								int nx, int ny, std::vector<color4> &cols)
		: ToolViewInterface(t, f,l2,h2,nx,ny), valy(vy) {

	std::vector<color4>::iterator it = cols.begin();
	for (;it!=cols.end();it++) {
		colors.push_back(*it);
	}


}


template <class T>
double ToolViewColor<T>::Get(int x, int y) {
	return valy[x*ndgl+y];
}



template <class T>
void ToolViewColor<T>::Draw(ToolPlot *p) {

	glBegin(GL_LINES);

	int satcount = colors.size(); //optim.sats.size();

	for (int i=0;i<ndis;i++) {
		float x = p->MapFloatToX(i*(xmax-xmin)/(ndis-1.)+xmin);

		int s = 0;
		for (int j=0;j<satcount;j++) {
			s += (int) Get(i,j);
		}

		if (s) {
			float y0 = p->MapToY(0);
			int s0 = 0;
			for (int j=0;j<satcount;j++) {
				s0 += (int) Get(i,j);
				float y = p->MapToY(s0/(float)s);
				colors[j](); //((glSatellite*)optim.sats[j])->col();
				//glColor4f(j*.3,1-j*.3,.8,.5);
				glVertex3f(x,-y0,0.3);
				glVertex3f(x,-y,0.3);

				y0 = y;
			}
		}
	}

	glEnd();

}

template <class T>
void ToolViewColor<T>::Text(std::ofstream &of) {


	int satcount = colors.size();
	for (int i=0;i<ndis;i++) {

		of << setw(20) << i*(xmax-xmin)/(ndis-1.)+xmin;

		for (int j=0;j<satcount;j++) {
			of << setw(20) << Get(i,j);
		}

		of << endl;
	}
}


template <class T>
void ToolViewColor<T>::GetMinMax(double &Low, double &High) {

	Low = 0;
	High = 1;

}



/*  --------------------------------------------- */



template class ToolViewLines<double>
;
template class ToolViewLines<int>
;
template class ToolViewBars<double>
;
template class ToolViewBars<int>
;
template class ToolViewMinMax<double,181>
;
template class ToolViewMinMax<int,181>
;
template class ToolViewColor<int>
;
