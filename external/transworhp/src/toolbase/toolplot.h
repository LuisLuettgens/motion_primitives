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
#ifndef TOOLPLOT_H
#define TOOLPLOT_H

#include "tool.h"

#include <vector>

#include "toolwindow.h"
#include "../base/minmax.h"

class ToolItem;
class ToolViewInterface;

class ToolPlot : public ToolWindow {
public:
	ToolPlot(const std::string &t, int w, int h, int &items);
	virtual ~ToolPlot();

	// void SetData(int l2,int h2, int *vy, int nx, int ny, int mod, int norm=0) ;
	void calcYscale();
	void Draw(bool isActive, int HT, const Point<int> &mouse, bool border=true);
	void Update() {}
	void Init() {}
	void Apply();
	void Reset();
	void Discard();
	virtual std::string GetFilename();
	void WriteText();
	void WriteImage();

	void drawFrame() const;
	void drawFrameText() const;
	void newHighLow();
	void DrawString(int x, int y, const std::string &s) const;
	//bool rawHighLow();

	float MapToY(double d) const;
	float MapFloatToX(double d) const;

	//protected:
	int LBorder,RBorder,TBorder,BBorder;
	double High,Low;
	//int High2,Low2;
	double yscale;
	void Finish();


	std::vector<ToolViewInterface*> views;
	std::vector<ToolViewInterface*>::iterator viewit;
	ToolViewInterface* currentView;

	std::vector<ToolItem*> items;
	int &itemsel;

void Resize();
double getStep() const;

protected:
	virtual bool mousebutton(int mx, int my, int button, int type) ;

};


#endif
