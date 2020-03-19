//
// C++ Implementation: toolbox
//
// Description:
//
//
// Author: Matthias Knauer,,, <tulio@visurgis>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifdef WIN32
#include "windows.h"
#endif

#include "GL/gl.h"

#include "conversion.h"
#include "tooltextbutton.h"


template <class T>
ToolTextbutton<T>::ToolTextbutton(const ToolRect &rect_, T &r, int id_)
		: ToolLabel<T>(rect_,r), id(id_) {}


template <class T>
void ToolTextbutton<T>::MouseClick(int button, int state, const Point<int> &p) {

	if (button==1 && state==SDL_MOUSEBUTTONDOWN) {
//std::cout << "WEWERW R" << std::endl;
		Message(id);
	
	}
}


template <class T>
void ToolTextbutton<T>::Draw(bool isActive, const Point<int> &mouse) {

int hilite = this->rect.contains(mouse);


	this->rect.border = Tool::OUTSET + ((hilite&1)?Tool::FILL : Tool::FILL2) + Tool::CENTER;
	ToolLabel<T>::Update();
	this->rect.Draw();
}


template class ToolTextbutton<const std::string>;

template class ToolTextbutton<std::string>;
