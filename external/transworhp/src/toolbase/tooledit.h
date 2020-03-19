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
#ifndef TOOLEDIT_H
#define TOOLEDIT_H

#include "tool.h"

#include "toolitem.h"



/** @ingroup toolbase
 *  @brief Edit-Elemente, z.B. Texte, Zahlen 
 */
template <class T>
class ToolEdit : public ToolData<T> {
public:
	ToolEdit(const ToolRect &rect_, T &r, int ai = 0, int *s=0)
			: ToolData<T>(rect_,r),applyimmediately(ai), scale(s) {
		edittext = GetString();
		cursor=edittext.size();
	}

	void Draw(bool isActive, const Point<int> &mouse);
	void MouseClick(int button, int state, const Point<int> &p);

	// Specialized...
	std::string GetString();
	void SetString(const std::string &s);

	void Key(SDL_Keysym &keysym);
	T ToValueType(const std::string& i) const;
	bool CheckKey(SDL_Keycode &c);

	bool TakesKey() {return true;}

	int cursor;
	std::string edittext; // Puffer zum Editieren
int applyimmediately;
	
	void Discard();
	
	int* scale;

};




#endif
