//
// C++ Interface: toolbox
//
// Description:
//
//
// Author: Matthias Knauer,,, <tulio@visurgis>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef TOOLVIEW_H
#define TOOLVIEW_H

#include "tool.h"

#include <vector>

#include "../base/minmax.h"

class ToolPlot;

class ToolViewInterface {
public:
	ToolViewInterface(const std::string &t, const std::string &f, int l2,int h2,
	                  int nx, int ny);
	virtual ~ToolViewInterface() {}

	virtual double Get(int x, int y) {
		return 0.;
	}
	virtual void Draw(ToolPlot *p) = 0;
	virtual void Text(std::ofstream &of) = 0;
	
	virtual void GetMinMax(double &Low, double &High);

	int xmin, xmax;
	int ndis, ndgl;
	

	std::string title;
	std::string filename;
	
};


class ToolViewGraph : public ToolViewInterface {
public:
	ToolViewGraph(const std::string &t, const std::string &f, int l2,int h2,std::vector<Point<double> >  *vy, const std::stringstream &str);

	void Draw(ToolPlot *p);
	void Text(std::ofstream &of);
	void GetMinMax(double &Low, double &High);

	std::vector<Point<double> > *valy;
std::string text;
};



template <class T, int N>
class ToolViewMinMax : public ToolViewInterface {
public:
	ToolViewMinMax<T,N>(const std::string &t, const std::string &f, int l2,int h2,MinMax<T,N>  *vy);

	double Get(int x, int y);
	void Draw(ToolPlot *p);
	void Text(std::ofstream &of);
	//void GetMinMax(double &Low, double &High);

	MinMax<T,N> *valy;
};


template <class T>
class ToolViewColor : public ToolViewInterface {
public:
	ToolViewColor<T>(const std::string &t, const std::string &f, int l2,int h2, T *vy,
	            int nx, int ny, std::vector<color4> &cols);

	double Get(int x, int y);
	void Draw(ToolPlot *p);
	void Text(std::ofstream &of);
	void GetMinMax(double &Low, double &High);

	T *valy;
	std::vector<color4> colors;

};



template <class T>
class ToolViewLines : public ToolViewInterface {
public:
	ToolViewLines<T>(const std::string &t, const std::string &f, int l2,int h2, T *vy, int nx);

	double Get(int x, int y);
	void Draw(ToolPlot *p);
	void Text(std::ofstream &of);
	//void GetMinMax(double &Low, double &High);

	T *valy;
};


template <class T>
class ToolViewBars : public ToolViewInterface {
public:
	ToolViewBars<T>(const std::string &t, const std::string &f, int l2,int h2, T *vy,
	            int nx, int nor);

	double Get(int x, int y);
	void Draw(ToolPlot *p);
	void Text(std::ofstream &of);
	void GetMinMax(double &Low, double &High);

	T *valy;
int norm;
};


#endif
