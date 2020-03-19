
#include <GL/glew.h>
#include <GL/gl.h>

#include "adjPlot.h"

#include <algorithm>

namespace tw {

/*
Matthias Rick
Klasse zum Ploten der Mu
TODO Hermite-Simpson
*/
adjPlot::adjPlot(int n, const double *zeit_, const double *mu_, int n_dis_, int n_ode_, const TWdiscretization *twdiscretization_) : BasePlot(0), zeit(zeit_), mu(mu_), n_dis(n_dis_), n_ode(n_ode_), komp(n), twdiscretization(twdiscretization_), plotData(2*(n_dis-1)) {

	//y Achse
	High = mu[komp];
	Low = mu[komp];

	//x Achse
	Low2 = zeit[0];
	High2 = zeit[n_dis - 1];
}

void adjPlot::drawData(DataStorage &/*ds*/, DataStorage &/*dstop*/, int /*ii*/) {

	if (mu) {

		// Werte zwischenspeichern
		for (int k = 0; k < n_dis-1; k++) {
			plotData[2*k] = static_cast<GLfloat>((zeit[k+1]+zeit[k])/2);
			plotData[2*k+1] = static_cast<GLfloat>(mu[k*n_ode*twdiscretization->punkte()+komp]);

			// Max. / Min. setzen
			if (plotData[2*k+1] < Low) Low = plotData[2*k+1];
			if (plotData[2*k+1] > High) High = plotData[2*k+1];
		}

		if (Low < -10000) {
			Low=-10000;
		}
		if (High > 10000) {
			High=10000;
		}

		const float farbe =  0.7f ;
		glColor3f(farbe + 0.5f, farbe, farbe);

		glPushMatrix();

		if (!IsIcon()) {
			glTranslatef(LBorder, TBorder, 0.0f);
			glScalef((GetWidth()-LBorder-RBorder)/(High2-Low2), scaledata*yscale, 1.0f);
			glTranslatef(0.0f, -Low, 0.0f);
		} else {
			glScalef(GetWidth()/(High2-Low2), scaledata*yscale, 1.0f);
		}

		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glEnableClientState(GL_VERTEX_ARRAY);

		glLineWidth(1.5f);

		glVertexPointer(2, GL_FLOAT, 0, plotData.data());
		glDrawArrays(GL_LINE_STRIP, 0,  n_dis-1);

		glLineWidth(1.0f);

		glDisableClientState(GL_VERTEX_ARRAY);

		glPopMatrix();
	}
}

void adjPlot::epsData(DataStorage &/*ds*/, DataStorage &/*dstop*/, EpsWriter */*epsw*/) {}

int adjPlot::rawHighLow(DataStorage &/*ds*/) {

	High = -1e20;
	Low = 1e20;

	if (mu) {
		for (int k = 0; k < n_dis-1; k++) {
			// Max. / Min. setzen
			double aux = mu[k*n_ode*twdiscretization->punkte()+komp];
			if (aux < Low) Low = aux;
			if (aux > High) High = aux;
		}
	}

	return 0;
}

}
