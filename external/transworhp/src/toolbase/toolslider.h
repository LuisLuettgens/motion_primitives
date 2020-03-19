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
#ifndef TOOLSLIDER_H
#define TOOLSLIDER_H

#include "tool.h"

#include "toolrect.h"
#include "toolitem.h"


/** @ingroup toolbase
 *  @brief Slider-Elemente, z.B. Texte, Zahlen
 */
template <class T>
class ToolSlider : public ToolData<T> {
public:
	ToolSlider(const ToolRect &rect_, T &r, T l, T h, int prec=4)
			: ToolData<T>(rect_,r), low(l), high(h), precision(prec) {
	}

	void Draw(bool isActive, const Point<int> &mouse);
	void MouseClick(int button, int state, const Point<int> &p);

	void Apply() {
	//	ref = value;
		//std::cout << "Item: " << ToString(ref) << std::endl;
	}
	
	// Specialized...
	std::string GetString();
	std::string GetRefString();
	void SetString(const std::string &s);

	void Key(SDL_Keysym &keysym);
	T ToValueType(const std::string& i) const;
	bool CheckKey(SDL_Keycode &c);

	bool TakesKey() {return true;}

	T low,high;
	
	int precision;
};


#endif
