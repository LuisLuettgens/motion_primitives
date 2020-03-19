//
// C++ Interface: plot
//
// Description:
//
//
// Author: Matthias Knauer <knauer@math.uni-bremen.de>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//

#pragma once

#include "baseplot.h"

#include "../core/TWparameter.h"

namespace tw {

class PhasePlot : public BasePlot {
public:
	PhasePlot(char c1, int d1, char c2, int d2, int ind);
	~PhasePlot();

protected:
	Selector data1,data2;

	int rawHighLow(DataStorage &ds) override;
	void drawData(DataStorage &ds, DataStorage &dstop, int ii) override;
	void epsData(DataStorage &ds, DataStorage &dstop, EpsWriter *epsw) override;
};


class MatrixPlot : public BasePlot {
public:
	MatrixPlot(double *d, int dim1,int dim2,int ind);
	~MatrixPlot();

protected:
	int D1, D2;
	double *matrix;

	int rawHighLow(DataStorage &ds) override;
	void drawData(DataStorage &ds, DataStorage &dstop, int ii) override;
	void epsData(DataStorage &ds, DataStorage &dstop, EpsWriter *epsw) override;
};

class DataPlot : public BasePlot {
public:
	DataPlot(std::string &s, Funktionenzeiger2 func, int par, int ind);
	DataPlot(std::string &s, Funktionenzeiger2 func, int *par, int ind);
	~DataPlot();

	void SetMaxTime(DataStorage::TimeMode_e timemode, double time) override;

protected:

	int rawHighLow(DataStorage &ds) override;
	void drawData(DataStorage &ds,DataStorage &dstop, int ii) override;
	void epsData(DataStorage &ds, DataStorage &dstop, EpsWriter *epsw) override;

	Funktionenzeiger2 Func;
	std::string title;
	//	int param;
	std::vector<int> indices;
};

/**
* Klasse zum Plotten des Fehlers und der Schrittweite
* @author Matthias Rick
*/
class gitterPlot : public BasePlot {
public:
	gitterPlot(std::string &s, int i, const std::vector<std::vector<double> > &zeit, const std::vector<std::vector<double> > &fehler);
	~gitterPlot();

protected:
	std::string title;
	/** Stuetzstellen */
	const std::vector<std::vector<double> > &zeit;
	/** Fehler */
	const std::vector<std::vector<double> > &fehler;

	/* fuer Achsen-Beschriftung */
	std::string ticktext2(double step, double i) const;
	void drawFrameText() const override;
	int rawHighLow(DataStorage &ds) override;
	void drawData(DataStorage &ds, DataStorage &dstop, int ii) override;
	void epsData(DataStorage &ds, DataStorage &dstop, EpsWriter *epsw) override;
};

/**
* Klasse zum Plotten der Stuetzstellen
* @author Matthias Rick
*/
class punktPlot : public BasePlot {
public:
	/** Erstellt einen Plot der Stuetzstellen
	* @param s Titel
	* @param i
	* @param zeit Diskretisierung
	* @param twdiscretization Diskretisierungstyp
	*/
	punktPlot(std::string &s, int i, const std::vector<std::vector<double> > &zeit, const TWdiscretization *twdiscretization);
	~punktPlot();

protected:
	/** Titel */
	std::string title;
	/** Stuetzstellen */
	const std::vector<std::vector<double> > &zeit;
	/** Diskretisierungstyp */
	const TWdiscretization *twdiscretization;

	int rawHighLow(DataStorage &ds) override;
	void drawData(DataStorage &ds, DataStorage &dstop, int ii) override;
	void epsData(DataStorage &ds, DataStorage &dstop, EpsWriter *epsw) override;
};

/**
* Klasse zum Plotten der Beschraenkungen
* @author Matthias Rick
*/
class lambdaPlot : public BasePlot {
public:
	/** Erstellt einen Plot der Beschraenkungen (Lambda)
	* @param s Titel
	* @param zeit Diskretisierung
	* @param lambda Multiplikatoren
	* @param n_dis Anzahl Stuetzstellen
	* @param n_ctrl Anzahl Steuerungen
	* @param n_ode Anzahl Zustaende
	* @param twdiscretization Diskretisierungstyp
	*/
	lambdaPlot(std::string &s, const double *zeit, const double *lambda, int n_dis, int n_ctrl, int n_ode, const TWdiscretization *twdiscretization);
	~lambdaPlot();

protected:
	/** Titel */
	std::string title;
	/** Stuetzstellen */
	const double *zeit;
	/** Multiplikatoren */
	const double *lambda;
	/** Anzahl Stutzstellen */
	int n_dis;
	/** Anzahl Steuerungen */
	int n_ctrl;
	/** Anzahl Zustaende */
	int n_ode;
	/** Diskretisierungstyp */
	const TWdiscretization *twdiscretization;

	int rawHighLow(DataStorage &ds) override;
	void drawData(DataStorage &ds, DataStorage &dstop, int ii) override;
	void epsData(DataStorage &ds, DataStorage &dstop, EpsWriter *epsw) override;
};

/**
* Klasse zum Plotten einer bestimmten Beschraenkung
* @author Matthias Rick
*/
class lambdaPlot2 : public lambdaPlot {

public:
	/** Erstellt einen Plot der Beschraenkungen (Lambda) fuer eine bestimme Komponente
	* @param s Titel
	* @param n Komponente
	* @param type Zustand=0 oder Steuerung=1
	* @param zeit Diskretisierung
	* @param lambda Multiplikatoren
	* @param n_dis Anzahl Stuetzstellen
	* @param n_ctrl Anzahl Steuerungen
	* @param n_ode Anzahl Zustaende
	* @param twdiscretization Diskretisierungstyp
	*/
	lambdaPlot2(std::string &s, int n, int type, const double *zeit, const double *lambda, int n_dis, int n_ctrl, int n_ode, const TWdiscretization *twdiscretization);
	~lambdaPlot2();

protected:
	/** Komponente (zB Zustand 3) */
	int komp;
	/** Typ (Zustand=0 oder Steuerung=1) */
	int type;
	/** Vector zum zwischenspeichern der Plot-Daten */
	std::vector<double> plotData;

	void drawData(DataStorage &ds, DataStorage &dstop, int ii) override;
};


class Data2Plot : public BasePlot {
public:
	Data2Plot(std::string& s, Funktionenzeiger2 func, int par1, int par2, int ind);
	~Data2Plot();

protected:

	int rawHighLow(DataStorage &ds) override;
	void drawData(DataStorage &ds, DataStorage &dstop, int ii) override;
	void epsData(DataStorage &ds, DataStorage &dstop, EpsWriter *epsw) override;

	Funktionenzeiger2 Func;
	std::string title;
	//	int param;
	std::vector<int> indices;
};

class UserPlot : public BasePlot {
public:
	UserPlot(Funktionenzeiger func,  int ind);
	UserPlot(FunktionenzeigerI func, int ind, int index_);
	UserPlot(FunktionenzeigerU func, int ind);
	~UserPlot();
	/*	int GetType() const {
			return 2;
		}
	*/
	void SetMaxTime2(double time) override;


	void SetUserControl(FunktionenzeigerC c) override;

	bool MouseInput(DataStorage &ds, int button, const Point<int> &p) override;

	void CallUserControl(DataStorage &ds, int button, double x, double y);

protected:
	Funktionenzeiger Func;
	FunktionenzeigerI FuncI;
	FunktionenzeigerC FuncC;
	FunktionenzeigerU FuncU;
	double maxtime;
	int lasttime;

	double userBorder[4];

	int rawHighLow(DataStorage &ds) override;
	void drawData(DataStorage &ds, DataStorage &dstop, int ii) override;
	void epsData(DataStorage &ds, DataStorage &dstop, EpsWriter *epsw) override;

	int iindex;
};

class TabularPlot : public BasePlot {
public:
	TabularPlot(int n, int m, double *val, char *heads, int ind);
	~TabularPlot();
	//	virtual void SetMaxTime(double time);

protected:
	int nn, mm;
	double *values;
	char *headers;

	std::string title;
	//	int param;

	int rawHighLow(DataStorage &ds) override;
	void drawData(DataStorage &ds, DataStorage &dstop, int ii) override;

	void epsData(DataStorage &ds, DataStorage &dstop, EpsWriter *epsw) override;

	void drawFrame() const override;
	void drawIconFrame() const override;
	void drawFrameText() const override;
	void drawTextData() const;
	void drawIconTextData() const;
};

}
