#pragma once

#include "baseplot.h"

namespace tw {

/** Klasse zum Plotten von Zustaenden und Steuerungen */
class GeneralPlot : public BasePlot {
public:
	GeneralPlot(char c, int d, int ind);
	~GeneralPlot();

	void drawData(DataStorage &ds, DataStorage &dstop, int ii) override;
	void SetMaxTime(DataStorage::TimeMode_e timemode, double time) override;
	void SetMinTime(DataStorage::TimeMode_e timemode, double time) override;

	int GetDgl() const;
	bool MouseInput(DataStorage &ds, int button,  const Point<int> &p) override;
	double MapFromX(DataStorage &ds, short x) const;

private:
	void epsData(DataStorage &ds, DataStorage &dstop, EpsWriter *epsw) override;
	int rawHighLow(DataStorage &ds) override;

	Selector data;
};

}
