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
#ifndef TOOLSCALE_H
#define TOOLSCALE_H

#include "tool.h"

#include "toolrect.h"
#include "toolitem.h"


/** @ingroup toolbase
 *  @brief Skala-Element
 */
template <class T>
class ToolScale : public ToolData<T> {
public:
	ToolScale(const ToolRect &rect_, T &r, double &min_, double &max_, int &step_, int&unitscale_)
			: ToolData<T>(rect_,r), minval(min_), maxval(max_), step(step_), unitscale(unitscale_) {

	}

	void Draw(bool isActive, const Point<int> &mouse);
	void MouseClick(int button, int state, const Point<int> &p);

	// Specialized...
	std::string GetString();
	void SetString(const std::string &s);

	void Key(SDL_Keysym &keysym);
	T ToValueType(const std::string& i) const;
	//bool CheckKey(SDLKey c);

	bool TakesKey() {return true;}

	double &minval,&maxval;
	int &step;
	int &unitscale;

private:
int getht(double i);
};

#endif
