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

using namespace std;

namespace tw {

/*
Matthias Rick
Klasse zum Ploten der Stuetzstellen
*/
punktPlot::punktPlot(string &s, int i, const vector<vector<double> > &zeit_, const TWdiscretization *twdiscretization_) : BasePlot(i), title(s), zeit(zeit_), twdiscretization(twdiscretization_) {
	//y Achse
	High = zeit.size();
	Low = 0;
	//x Achse
	High2 = 1;
	Low2 = 0;
}

punktPlot::~punktPlot() {}

/*
Methode zum zeichnen
*/
void punktPlot::drawData(DataStorage &/*ds*/, DataStorage &/*dstop*/, int /*ii*/) {
	
	float farbe = 0.0f;

	glColor3f(farbe, farbe, farbe);
	glPointSize(1);

	glBegin(GL_POINTS);

	for (size_t j = 0; j < zeit.size(); j++) {
		
		float y1 = MapToY(j);
		
		for (size_t k = 0; k < zeit[j].size(); k++) {
			float x1 = MapFloatToX(zeit[j][k]);
			glVertex3d(x1, y1, 0.5f);
			
			// Hermite-Simpson: Mittelpunkte plotten
			if (twdiscretization->type == TWdiscretizationType::HermiteSimpson && k != zeit[j].size()-1) {
				glColor3f(0.5f,0.5f,0.5f);
				x1 = MapFloatToX((zeit[j][k]+zeit[j][k+1])*0.5);
				glVertex3d(x1, y1, 0.5f);
				glColor3f(farbe, farbe, farbe);
			}
		}
	}

	glEnd();
	
	//damit die folgenden Plots nicht veraendert werden
	glPointSize(1);
}

void punktPlot::epsData(DataStorage &/*ds*/, DataStorage &/*dstop*/, EpsWriter */*epsw*/) {
}

int punktPlot::rawHighLow(DataStorage &/*ds*/) {
	return 0;
}

}
