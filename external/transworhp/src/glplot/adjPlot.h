#pragma once

#include "baseplot.h"

#include "../core/TWparameter.h"

namespace tw {

/**
* Klasse zum Plotten der Adjungierten (Mu)
* @author Matthias Rick
*/
class adjPlot : public BasePlot {

public:
	/** Erstellt einen Plot der Adjungierten (Mu)
	* @param n Komponente
	* @param zeit Diskretisierung
	* @param mu Adjungierten
	* @param n_dis Anzahl Stuetzstellen
	* @param n_ctrl Anzahl Steuerungen
	* @param n_ode Anzahl Zustaende
	* @param twdiscretization Diskretisierungstyp
	*/
	adjPlot(int n, const double *zeit, const double *mu, int n_dis, int n_ode, const TWdiscretization *twdiscretization);

private:
	/** Stuetzstellen */
	const double *zeit;
	/** Adjungierten */
	const double *mu;
	/** Anzhal Stuetzstellen */
	int n_dis;
	/** Anzahl Zustaende */
	int n_ode;
	/** Komponente */
	int komp;
	/** Diskretisierungstyp */
	const TWdiscretization *twdiscretization;
	/** Vector zum zwischenspeichern der Plot-Daten */
	std::vector<GLfloat> plotData;

	void drawData(DataStorage &ds, DataStorage &dstop, int ii) override;
	void epsData(DataStorage &ds, DataStorage &dstop, EpsWriter *epsw) override;
	int rawHighLow(DataStorage &ds) override;
};

}
