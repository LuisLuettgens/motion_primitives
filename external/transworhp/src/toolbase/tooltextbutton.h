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
#ifndef TOOLTEXTBUTTON_H
#define TOOLTEXTBUTTON_H

#include "tool.h"

#include "toolrect.h"
#include "toolitem.h"


/** @ingroup toolbase
 *  @brief ToolTextbutton.
 *
 */
template <class T>
class ToolTextbutton : public ToolLabel<T> {
public:
	ToolTextbutton(const ToolRect &rect_, T &r, int id_);

	void Draw(bool isActive, const Point<int> &mouse);
	void MouseClick(int button, int state, const Point<int> &p);

	bool TakesKey() {return true;}
	
private:
	int id;

};

#endif
