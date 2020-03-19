//
// Author: Matthias Knauer <knauer@math.uni-bremen.de>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
#pragma once

#include "../base/point.h"
#include "SDL2/SDL.h"

#include <vector>

namespace tw {

class ToolMenu;

/** @ingroup toolbase
 *  @brief MenuEntry
 *
 */
class ToolMenuEntry {

public:
	ToolMenuEntry() : altkey(SDLK_CLEAR), id(0), openentry2(0) {}
	ToolMenuEntry(const std::string &t, int i, const std::string &hot, int *sel);
	virtual ~ToolMenuEntry();

	virtual void drawText(int x, int y, float z, int ww, int totalwid, int open,Point<int> &p1,Point<int> &p2);
	void DrawOpenMenu(const ToolMenu *m, int x, int y, float z);
	bool MouseMotion(const Point<int> &m,int x,int y);
	bool MouseButton(const Point<int> &mouse, int x, int y, Uint8 type, Uint8 button);
	virtual bool MouseButton1(const Point<int> &mouse, int x, int y, Uint8 type, Uint8 button);
	int GetMenuWidth();
	int GetMenuHeight();
	virtual bool KeyBoardFunc(SDL_Keycode key);
	void ResetSelection();

public:
	std::string text;
	SDL_Keycode altkey;

	std::string hottext;
	int id;

	std::vector<ToolMenuEntry*> me;
	int pos;
	int w,w2,w3;
	int linefrom,lineto;

	ToolMenuEntry *openentry2;
	int yy;

	int *selected;

};

class ToolMenuSliderEntry : public ToolMenuEntry {
public:
	ToolMenuSliderEntry(const std::string &t, int i, const std::string &hot, float *val, float minval, float maxval);

	void drawText(int x, int y, float z, int ww, int totalwid, int open,Point<int> &p1,Point<int> &p2);
	bool KeyBoardFunc(SDL_Keycode key);
	bool MouseButton1(const Point<int> &mouse, int x, int y, Uint8 type, Uint8 button);
	float *val;
	float minval, maxval;
};

}
