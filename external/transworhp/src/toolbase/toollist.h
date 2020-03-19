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
#ifndef TOOLLIST_H
#define TOOLLIST_H

#include "tool.h"

#include "toolrect.h"
#include "toolitem.h"

/** @ingroup toolbase
 *  @brief Listenelemente, z.B. Text-Aufzählungen oder Farben 
 */
 

template <class T>
class ToolList : public ToolData<T> {

public:
/** Constructor */
	ToolList(const ToolRect &rect_, T &r, std::string *Map);
	ToolList(const ToolRect &rect_, T &r, const std::vector<std::string> &lst);
	virtual ~ToolList();
	
	void Open();
	void Close();

	void MouseClick(int button, int state, const Point<int> &p);
	void Draw(bool isActive, const Point<int> &mouse);
	void Key(SDL_Keysym &keysym);

	// Specialized for int and Color4
	std::string GetString();
	void SetString(const std::string &s);

	void SetValue(int t);
	int GetValue();
	T ToValueType(int i) const;

	bool TakesKey() {return true;}

private:
	void drawItem(ToolRect &r, T i, double b, float addz=0);

	std::string *map;
	int opened;
	ToolRect openRect;
	ToolRect closeRect;
	int size;
	int index;
	int deleting;
	
};

#endif
