#ifdef WIN32
#include "windows.h"
#endif
#include "plot.h"
//#include "gldisplay.h"
#include <GL/gl.h>

#include "xopt_eps.h"

#include <iostream>
#include <iomanip>

#include "../toolbase/tool.h"

#ifdef WIN32
#define fabs(A) ((A)>0?(A):-(A))
#ifdef MINGW
#define min(A,B) ((A)<(B)?(A):(B))
#define max(A,B) ((A)>(B)?(A):(B))
#endif
#else
#define min(A,B) ((A)<(B)?(A):(B))
#define max(A,B) ((A)>(B)?(A):(B))
#define fabs(A) ((A)>0?(A):-(A))
#endif

using namespace std;

namespace tw {

/*
Matthias Rick
Klasse zum Ploten der Schrittweite und des Diskretisierungsfehlers
i == 0 Zeit (Schrittweite)
i == 1 Fehler
*/
gitterPlot::gitterPlot(string &s, int i, const vector<vector<double> > &zeit_, const vector<vector<double> > &fehler_) : BasePlot(i), title(s), zeit(zeit_), fehler(fehler_) {
	//y Achse
	//High = 2;
	High = 0;
	Low = 0;
	if (index == 0) {
		//Low = -8;
		Low = -5;
	} else {
		Low = -10;
		//Low = -16;
	}

	//x Achse
	High2 = 1;
	Low2 = 0;
}

gitterPlot::~gitterPlot() {}

/*
Methode zum zeichnen
*/
void gitterPlot::drawData(DataStorage &/*ds*/, DataStorage &/*dstop*/, int /*ii*/) {
	double x1, x2, y1, y2;
	double farbe;

	y1 = MapToY(Low);

	if (index == 0) { // Schrittweite

		if (zeit.size() == 1) {
			High = (int) log10(zeit[0][1] - zeit[0][0]);
		}

		for (size_t j = 0; j < zeit.size(); j++) {


			farbe = ((double) 1 / zeit.size()) * j;
			glColor3f(farbe + 0.5, farbe, farbe);

			for (size_t k = 1; k < zeit[j].size(); k++) {

				//cout << "j:" << j << ", k:" << k << ", " << zeit[j][k] << endl;

				y2 = MapToY(log10(zeit[j][k] - zeit[j][k - 1]));

				x1 = MapFloatToX(zeit[j][k]);
				x2 = MapFloatToX(zeit[j][k - 1]);

				// Umrandung
				if (j == zeit.size() - 1) {
					glColor3f(0, 0, 0);

					glBegin(GL_TRIANGLE_FAN);

					glVertex3d(x1 + 1, y1 - 1, j - 0.5);
					glVertex3d(x1 + 1, y2 + 1, j - 0.5);
					glVertex3d(x2 - 1, y2 + 1, j - 0.5);
					glVertex3d(x2 - 1, y1 - 1, j - 0.5);

					glEnd();

					glColor3f(farbe+0.5, farbe, farbe);
				}


				glBegin(GL_TRIANGLE_FAN);

				glVertex3d(x1, y1, j);
				glVertex3d(x1, y2, j);
				glVertex3d(x2, y2, j);
				glVertex3d(x2, y1, j);

				glEnd();
			}
		}
	}
	else if (index == 1) { // Fehler

		for (size_t j = 0; j < fehler.size(); j++) {

			farbe = ((double) 1 / zeit.size()) * j;
			glColor3f(farbe, farbe, farbe);

			for (size_t k = 1; k < fehler[j].size()+1; k++) {

				//cout << "j:" << j << ", k:" << k << ", " << fehler[j][k] << endl;

				if (log10(fehler[j][k - 1]) >= -10) {
					y2 = MapToY(log10(fehler[j][k - 1]));
				}
				else {
					y2 = y1;
				}

				//if (k >= zeit[j].size())
				//	cout << "++++++++********++++++------>>>> " << k << " " << zeit[j].size() << endl;

				x1 = MapFloatToX(zeit[j][k]);
				x2 = MapFloatToX(zeit[j][k - 1]);

				// Umrandung
				if (j == fehler.size() - 1) {
					glColor3f(0, 0, 0);

					glBegin(GL_TRIANGLE_FAN);

					glVertex3d(x1+1, y1-1, j - 0.5);
					glVertex3d(x1+1, y2+1, j - 0.5);
					glVertex3d(x2-1, y2+1, j - 0.5);
					glVertex3d(x2-1, y1-1, j - 0.5);

					glEnd();

					glColor3f(farbe, farbe, farbe);
				}

				glBegin(GL_TRIANGLE_FAN);

				glVertex3d(x1, y1, j);
				glVertex3d(x1, y2, j);
				glVertex3d(x2, y2, j);
				glVertex3d(x2, y1, j);

				glEnd();
			}
		}

	}

}

std::string gitterPlot::ticktext2(double step, double i) const {
	std::stringstream s;

	if (step>=10000) {
		//s.setf(ios::left, ios::adjustfield);
		s.setf(std::ios::scientific, ios::floatfield);
		int ii = i;
		s << std::setprecision(2) << std::setw(10) << "10^" << (double)ii;
	} else if (step>=1)
		s << std::setw(5) << "10^" << (int)i;
	else if (step<.0001)
		s << std::setiosflags(std::ios::showpoint | std::ios::right |
		                      std::ios::fixed)
		  << std::setprecision(5) << std::setw(7) << "10^" << i;
	else if (step<.001)
		s << std::setiosflags(std::ios::showpoint | std::ios::right |
		                      std::ios::fixed)
		  << std::setprecision(4) << std::setw(6) << "10^" << i;
	else if (step<.01)
		s << std::setiosflags(std::ios::showpoint | std::ios::right |
		                      std::ios::fixed)
		  << std::setprecision(3) << std::setw(5) << "10^" << i;
	else if (step<.1)
		s << std::setiosflags(std::ios::showpoint | std::ios::right |
		                      std::ios::fixed)
		  << std::setprecision(2) << std::setw(5) << "10^" << i;
	else
		s << setiosflags(std::ios::showpoint | std::ios::right | std::ios::fixed)
		 << std::setprecision(1) << std::setw(5) << "10^" << i;

	return s.str();
}

void gitterPlot::drawFrameText() const {

	int width = GetWidth();
	int height = GetHeight();


	double step = 25*(High2-Low2)/(width-RBorder-LBorder);

	step = roundStep(step);


	for (double i=0; i<=High2+.001; i+=step) {
		float xzero=MapFloatToX(i);


		if (xzero>=LBorder-2 && xzero<=width-RBorder+2) {
			std::string buf = ticktext(step,i);
			//			SetForeground(grey);
			//			if (buf=="2")
			//				int aaa = 0;

			int ww = Tool::font->StringLength(buf.c_str());

			//int pos = xzero-buf.length()*3+1;
			int pos = xzero-ww/2;
			DrawString(pos,TBorder-10-3,buf);
		}
	}

	for (double i=-step; i>=Low2-.001; i-=step) {
		float xzero=MapFloatToX(i);

		if (xzero>=LBorder-2 && xzero<=width-RBorder+2) {
			std::string buf = ticktext(step,i);
			//	SetForeground(grey);
			//int pos = xzero-buf.length()*3+1;

			int ww = Tool::font->StringLength(buf.c_str());
			int pos = xzero-ww/2;
			DrawString(pos,TBorder-10-3, buf);
		}
	}


	step = 15*(High-Low)/(height-TBorder-BBorder);
	step = roundStep(step);

	double l = max(0., Low);
	l = (int)(l/step)*step;

	for (double i=l; i<=High+.001; i+=step) {
		float yzero=MapToY(i);

		if (yzero>=TBorder-2 && yzero<=height-BBorder+2) {
			//		SetForeground(grey);

			std::string buf = ticktext2(step,i);

			int ww = Tool::font->StringLength(buf.c_str());
			int pos = LBorder-ww-2;

			DrawString(pos,yzero-5,buf);
		}
	}

	l = min(0., High);
	l = (int)(l/step)*step;
	if (l>=0.)
		l= -step;


	for (double i=l; i>=Low-.001; i-=step) {
		float yzero=MapToY(i);

		if (yzero>=TBorder-2 && yzero<=height-BBorder+2) {

			std::string buf = ticktext2(step,i);

			int ww = Tool::font->StringLength(buf.c_str());
			int pos = LBorder-ww-2;
			//		SetForeground(grey);

			DrawString(pos,yzero-5,buf);
		}
	}
}

void gitterPlot::epsData(DataStorage &/*ds*/, DataStorage &/*dstop*/, EpsWriter */*epsw*/) {
}

int gitterPlot::rawHighLow(DataStorage &/*ds*/) {
	return 0;
}

}
