#ifndef TOOLFLOATING_H
#define TOOLFLOATING_H

#include "toolrect.h"

class ToolItem;

class ToolFloating {

public:
	ToolFloating();
	~ToolFloating();
	bool MouseMotion(SDL_MouseMotionEvent &m);
	bool MouseButton(SDL_MouseButtonEvent &m);

	void Draw(int HT, int WD, const Point<int> &mouse);

	ToolItem *item;

private:
	int moving;

	Point<int> mouse;
	ToolRect rect0;

};


#endif
