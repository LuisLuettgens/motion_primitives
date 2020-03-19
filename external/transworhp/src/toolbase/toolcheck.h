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
#ifndef TOOLCHECK_H
#define TOOLCHECK_H

#include "tool.h"

#include "toolrect.h"
#include "toolitem.h"

/** @ingroup toolbase
 *  @brief Check-Elemente, z.B. Texte, Zahlen
 */
template <class T>
class ToolCheck : public ToolData<T> {
public:
	ToolCheck(const ToolRect &rect_, T &r, const std::string &t)
			: ToolData<T>(rect_,r) {

		this->rect.text = t;	
		
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


//	std::string text;

};




#endif
