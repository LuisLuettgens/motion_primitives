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
#ifndef TOOLBUTTON_H
#define TOOLBUTTON_H

#include "tool.h"
#include <vector>

#include "toolwindow.h"

#include "../base/dynloader.h"

namespace tw {

/** @ingroup toolbase
 *  @brief ToolButton.
 *
 */
class ToolButton {

public:
	/** Constructor. */
	ToolButton(int id_, int image_, int x_,int y_);

	/** Destructor. */
	virtual ~ToolButton() {}

	virtual void Draw(const Point<int> &mouse);
	bool contains(int x_,int y_);
	virtual float dim();
	virtual int imagesel();
	virtual int pressed() {return 0;}
	
	int id;
	int imagex,imagey;
	int x,y;
	int width,height;
	int col;

};


/** @ingroup toolbase
 *  @brief ToolStateButton.
 *
 */
class ToolStateButton : public ToolButton {
public:

	/** Constructor. */
	ToolStateButton(int &s, int id_, int image_, int x_,int y_, int modes_);

	float dim();
	int imagesel();
	int pressed();	
	int &state;
	int modes;
};

}

#endif
