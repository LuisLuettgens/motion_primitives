
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

#include "../core/TWparameter.h"
#include <algorithm>

using namespace std;

namespace tw {

/*
Matthias Rick
Klasse zum Ploten der Lambda
*/
lambdaPlot::lambdaPlot(std::string &s, const double *zeit_, const double *lambda_, int n_dis_, int n_ctrl_, int n_ode_, const TWdiscretization *twdiscretization_) : BasePlot(0), title(s), zeit(zeit_), lambda(lambda_), n_dis(n_dis_), n_ctrl(n_ctrl_), n_ode(n_ode_), twdiscretization(twdiscretization_) {
	//y Achse
	High = 1;
	Low = -1;

	//x Achse
	High2 = 1;
	Low2 = 0;
}

lambdaPlot::~lambdaPlot() {}

void lambdaPlot::drawData(DataStorage &/*ds*/, DataStorage &/*dstop*/, int /*ii*/) {
	double x1, y1;
	double farbe;

	// TODO n_param
	High = *std::max_element(lambda,lambda+(n_ode + n_ctrl) * twdiscretization->stuetzstellen(n_dis));// + n_param;
	                                     //(n_ode + n_ctrl) * twdiscretization->punkte() * n_dis + n_ode);
	
	Low = *std::min_element(lambda,lambda+(n_ode + n_ctrl) * twdiscretization->stuetzstellen(n_dis));// + n_param;
					    //(n_ode + n_ctrl) * twdiscretization->punkte() * n_dis + n_ode);
	if (Low<-10000) Low=-10000;
	if (High>10000) High=10000;

	glPointSize(3.);

	for (int k = 0; k < n_dis; k++) {

		glBegin(GL_POINTS);
		for (int ctrl = 0; ctrl < n_ctrl; ctrl++) {
		  
			farbe = ((double) 1.0 / n_ctrl) * ctrl;
			glColor3f(farbe + 0.8, farbe, farbe);

			if (abs(lambda[(n_ode + n_ctrl) * twdiscretization->punkte() * k + n_ode + ctrl]) > 0) {
				y1 = MapToY(lambda[(n_ode + n_ctrl) * twdiscretization->punkte() * k + n_ode + ctrl]);
				x1 = MapFloatToX(zeit[k]);
				glVertex3d(x1, y1, 0.5);
			}
		}
		glEnd();
		
		glBegin(GL_POINTS);
		for (int ode = 0; ode < n_ode; ode++) {
		  
			farbe = ((double) 1.0 / n_ode) * ode;
			glColor3f(farbe, farbe, farbe+0.8);

			if (abs(lambda[(n_ode + n_ctrl) * twdiscretization->punkte() * k + ode]) > 0) {
				y1 = MapToY(lambda[(n_ode + n_ctrl) * twdiscretization->punkte() * k + ode]);
				x1 = MapFloatToX(zeit[k]);
				glVertex3d(x1, y1, 0.5);
			}
		}
		glEnd();
	}
	
	glPointSize(1.);
	glLineWidth(0.1);
	
	glEnable(GL_LINE_SMOOTH);
	
	for (int ctrl = 0; ctrl < n_ctrl; ctrl++) {

		farbe = ((double) 1.0 / n_ctrl) * ctrl;
		glColor3f(farbe+0.8, farbe, farbe);
		
		glBegin(GL_LINE_STRIP);
		for (int k = 0; k < n_dis; k++) {
			y1 = MapToY(lambda[(n_ode + n_ctrl) * twdiscretization->punkte() * k + n_ode + ctrl]);
			x1 = MapFloatToX(zeit[k]);
			glVertex3f(x1, y1, 0.5);
		}
		glEnd();
		
		if (twdiscretization->type == TWdiscretizationType::HermiteSimpson) {
			glBegin(GL_LINE_STRIP);
			for (int k = 0; k < n_dis-1; k++) {
				y1 = MapToY(lambda[(n_ode + n_ctrl) * twdiscretization->punkte() * k + n_ode + (n_ode + n_ctrl) + ctrl]);
				x1 = MapFloatToX((zeit[k]+zeit[k+1])*.5);
				glVertex3f(x1, y1, 0.5);
			}
			glEnd();
		}
	}
	
	for (int ode = 0; ode < n_ode; ode++) {
	  
		farbe = ((double) 1.0 / n_ode) * ode;
		glColor3f(farbe, farbe, farbe+0.8);
	  
		glBegin(GL_LINE_STRIP);
		for (int k = 0; k < n_dis; k++) {
			y1 = MapToY(lambda[(n_ode + n_ctrl) * twdiscretization->punkte() * k + ode]);
			x1 = MapFloatToX(zeit[k]);
			glVertex3d(x1, y1, 0.5);
		}
		glEnd();
		
		if (twdiscretization->type == TWdiscretizationType::HermiteSimpson) {
			//glBegin(GL_LINE_STRIP);
			for (int k = 0; k < n_dis-1; k++) {
				y1 = MapToY(lambda[(n_ode + n_ctrl) * twdiscretization->punkte() * k + (n_ode + n_ctrl) + ode]);
				x1 = MapFloatToX((zeit[k]+zeit[k+1])*.5);
				glVertex3f(x1, y1, 0.5);
			}
			//glEnd;
		}
	}
	
	glDisable(GL_LINE_SMOOTH);
	
	glLineWidth(1.);
	
}

void lambdaPlot::epsData(DataStorage &/*ds*/, DataStorage &/*dstop*/, EpsWriter */*epsw*/) {
}

int lambdaPlot::rawHighLow(DataStorage &/*ds*/) {
	return 0;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// Klasse zum Plotten einzelner Komponente ///////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////


lambdaPlot2::lambdaPlot2(std::string &s, int n_, int type_, const double *zeit_, const double *lambda_, int n_dis_, int n_ctrl_, int n_ode_, const TWdiscretization *twdiscretization_) : lambdaPlot(s,zeit_,lambda_,n_dis_,n_ctrl_,n_ode_,twdiscretization_), komp(n_), type(type_) {
	
	plotData.reserve(n_dis);
	plotData.resize(n_dis);
}

lambdaPlot2::~lambdaPlot2() {}

//TODO Hermite-Simpson
void lambdaPlot2::drawData(DataStorage &/*ds*/, DataStorage &/*dstop*/, int /*ii*/) {
 

	for (int k = 0; k < n_dis; ++k) {
		if (type == 0) {
			plotData[k] = lambda[(n_ode + n_ctrl) * twdiscretization->punkte() * k + komp];
		} else if (type == 1) {
			plotData[k] = lambda[(n_ode + n_ctrl) * twdiscretization->punkte() * k + n_ode + komp];
		}
	}
	
	// Max und Min bestimmen
	High = *std::max_element(plotData.begin(),plotData.end());
	Low = *std::min_element(plotData.begin(),plotData.end());
	
	if (Low == High) {
		Low -= 0.1;
		High += 0.1;
	}
	
	if (Low<-10000) Low=-10000;
	if (High>10000) High=10000;

	glPointSize(3.);

	for (int k = 0; k < n_dis; k++) {

		if (type == 1) {
			glBegin(GL_POINTS);
			
			double farbe = ((double) 1.0 / n_ctrl) * komp;
			glColor3f(farbe + 0.8, farbe, farbe);

			if (abs(plotData[k]) > 0) {
				double y1 = MapToY(plotData[k]);
				double x1 = MapFloatToX(zeit[k]);
				glVertex3d(x1, y1, 0.5);
			}
			
			glEnd();
		}
		
		else if (type == 0) {
			glBegin(GL_POINTS);
			
			double farbe = ((double) 1.0 / n_ode) * komp;
			glColor3f(farbe, farbe, farbe+0.8);

			if (abs(plotData[k]) > 0) {
				double y1 = MapToY(plotData[k]);
				double x1 = MapFloatToX(zeit[k]);
				glVertex3d(x1, y1, 0.5);
			}
			
			glEnd();
		}
	}
	
	glPointSize(1.);
	glLineWidth(0.1);
	
	glEnable(GL_LINE_SMOOTH);
	
	if (type == 1) {

		double farbe = ((double) 1.0 / n_ctrl) * komp;
		glColor3f(farbe+0.8, farbe, farbe);
		
		glBegin(GL_LINE_STRIP);
		for (int k = 0; k < n_dis; k++) {
			double y1 = MapToY(plotData[k]);
			double x1 = MapFloatToX(zeit[k]);
			glVertex3f(x1, y1, 0.5);
		}
		glEnd();
		
		/*
		if (twdiscretization->type == TWdiscretizationType::HermiteSimpson) {
			glBegin(GL_LINE_STRIP);
			for (int k = 0; k < n_dis-1; k++) {
				double y1 = MapToY(lambda[(n_ode + n_ctrl) * twdiscretization->punkte() * k + n_ode + (n_ode + n_ctrl) + komp]);
				double x1 = MapFloatToX((zeit[k]+zeit[k+1])*.5);
				glVertex3f(x1, y1, 0.5);
			}
			glEnd();
		}
		*/
	}
	
	else if (type == 0) {
		
		double farbe = ((double) 1.0 / n_ode) * komp;
		glColor3f(farbe, farbe, farbe+0.8);
		
		glBegin(GL_LINE_STRIP);
		for (int k = 0; k < n_dis; k++) {
			double y1 = MapToY(plotData[k]);
			double x1 = MapFloatToX(zeit[k]);
			glVertex3d(x1, y1, 0.5);
		}
		glEnd();
		
		/*
		if (twdiscretization->type == TWdiscretizationType::HermiteSimpson) {
			glBegin(GL_LINE_STRIP);
			for (int k = 0; k < n_dis-1; k++) {
				double y1 = MapToY(lambda[(n_ode + n_ctrl) * twdiscretization->punkte() * k + (n_ode + n_ctrl) + komp]);
				double x1 = MapFloatToX((zeit[k]+zeit[k+1])*.5);
				glVertex3f(x1, y1, 0.5);
			}
			glEnd;
		}
		*/
	}
	
	glDisable(GL_LINE_SMOOTH);
	
	glLineWidth(1.);
	
}

}
