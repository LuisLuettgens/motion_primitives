#pragma once

#include "plot.h"

#include <string>
#include <fstream>
#include <vector>
#include <map>

namespace tw {

/** @defgroup eps Writing EPS files
 *  ???
 * @{
 */
class EpsWriter {
public:

	EpsWriter(const std::string &filename);
	~EpsWriter();

	void Header(const std::string &Title);
	void Footer();
	void Frame(const std::string &title, BasePlot *pl);

	void SetLineColor(TWcolor color);
	void SetTextColor(TWcolor color);
	void LineWidth(float width);

	void SetLineColorI(int  color);

	// einfache Linien und Texte
	void Dot(float x1, float y1, float sz=4.f);
	void Line(float x1, float y1, float x2, float y2);
	void TextAlign(float x, float y, const std::string &text,
	               int xalign, int yalign);

	// mehrfache Linien und Texte
	void Rectangle(float x1, float y1, float x2, float y2, bool clip=false);

	void XLine(const std::vector<float> &x, float y1, float y2);
	void YLine(const std::vector<float> &y, float x1, float x2);

	void XTextAlign(const std::map<float, std::string> &x, float y,
	                int xalign, int yalign);
	void YTextAlign(const std::map<float, std::string> &y, float x,
	                int xalign, int yalign);


	// Punkte in Array abspeichern und alles malen.
	void AddPoint(float x, float y);
	void Dot(float sz=4.f);
	void Line();
	void ClearPoints();


	float epsMapToY(double d) const;
	float epsMapToX(int x, int ndis) const;
	float epsMapFloatToX(double d) const;
	float epsMapTimeToX(DataStorage &ds, int i);
	float epsMapStopTimeToX(DataStorage &ds, DataStorage &dstop, int i) const;

	//protected:


	bool Mirror;

	int Width2;
	int Height2;

	int _LBorder;
	/** Right border */
	int _RBorder;
	/** Upper border */
	int _TBorder;
	/** Lower border */
	int _BBorder;


	bool absetzen;

	double High, Low, High2, Low2;



private:
	float linewidth;
	std::vector<float> xx;
	std::vector<float> yy;
	std::string textColor;
	std::string lineColor;

	std::ofstream file;

	};

/** @} */

extern EpsWriter *globEpswriter;

double roundHiLo(double diff);
std::string ticktext(double step, double i);
std::string ticktext2(double step, double i);
extern double stopmark;

}
