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
#include "toolbox.h"

#include "GL/gl.h"
#include "conversion.h"
#include "../base/point.h"
#include "toolbutton.h"

using namespace std;


ToolBox::ToolBox(const std::string &t) : ToolWindow(t) {

	AddIcon(ToolIcon::OPEN);

	hasactiveitem=false;

	
}


ToolBox::~ToolBox() {
	clearpointer(buttons);
	clearpointer(items);

}

void ToolBox::Draw(bool isActive, int HT, const Point<int> &mouse, bool border) {

	if (!visible)
		return;

	// Iconleiste

	glLineWidth(1);
	glDisable(GL_LINE_SMOOTH);
	
	ToolWindow::Draw(isActive, HT, mouse, border);


	// Projektion f�r inneres (mit Clip + Scroll
	innerProj(HT);


	glLineWidth(1);

	//int mx = m.x-x;
	//	int my = m.y-y-theight;

	//	mousemotion(mx,my+scrollcur);

	Point<int> m1 = mouse - Point<int>(x,y+theight+22-scrollcur);
	Point<int> m2 = mouse - Point<int>(x,y+theight-scrollcur);

	if (open) {
		std::vector<ToolButton*>::const_iterator it = buttons.begin();
		int i=0;
		for (;it!=buttons.end();it++,i++) {
			Tool::colors[0].Dim(moving?.3:.2)();
			(*it)->Draw(m1);
		}

		Update();


		std::vector<ToolItem*>::const_iterator it2 = items.begin();
		for (;it2!=items.end();it2++) {

			bool act=false;
			if (hasactiveitem)
				act = (activeItem == it2);

			(*it2)->Draw(act && isActive, m2);
		}


		Tool::colors[6].Alpha(.6)();
		drawItem();

	}
	glLineWidth(1);

}

void ToolBox::Apply() {

	std::vector<ToolItem*>::const_iterator it2 = items.begin();
	for (;it2!=items.end();it2++) {
		(*it2)->Apply();
	}

}

void ToolBox::Discard() {

	std::vector<ToolItem*>::const_iterator it2 = items.begin();
	for (;it2!=items.end();it2++) {
		(*it2)->Discard();
	}

}

void ToolBox::Init() {

	std::vector<ToolItem*>::const_iterator it2 = items.begin();
	for (;it2!=items.end();it2++) {
		(*it2)->Discard();
	}
	
	//activeItem = items.begin();
	hasactiveitem=false;
	activeItem = items.end();
}



/*bool ToolBox::MouseMotion(SDL_MouseMotionEvent &m) {
	if (moving) {
		int tmpx = m.x-mouse.x;
		int tmpy = m.y-mouse.y;
 
		x = tmpx;
		y = tmpy;
		return true;
	}
	return false;
}
*/

bool ToolBox::mousemotion(int mx, int my) {

	std::vector<ToolButton*>::iterator it = buttons.begin();
	for (;it!=buttons.end();it++) {

		if ((*it)->contains(mx,my-theight-4)) {
			//			hiliteButton = (*it);
			//	hiliteItem = 0;
			return false;
		}
	}

	//	cout << mx << " " << my << endl;

	std::vector<ToolItem*>::iterator it2 = items.begin();
	for (;it2!=items.end();it2++) {
		//if ((*it2)->TakesKey()) {
		if ((*it2)->contains(mx,my)) {
			//(*it2)->MouseClick(button,type,Point(mx,my));
			//		hiliteItem = (*it2);
			//	hiliteButton = 0;
			return false;
		}
	}

	//hiliteItem = 0;
	//hiliteButton = 0;
	return false;
}


bool ToolBox::mousebutton(int mx, int my, int button, int type) {

	/*bool ToolBox::MouseButton(SDL_MouseButtonEvent &m) {
	 
		if ((m.button == 1 || m.button==4 || m.button==5) && m.type == SDL_MOUSEBUTTONUP) {
			if (moving)
				moving=0;
		}
	 
		if ((m.button == 1 || m.button==4 || m.button==5) && m.type == SDL_MOUSEBUTTONDOWN) {
	 
	// in Fenster
			if (m.x>x && m.x<x+width && m.y>y && m.y<y+theight) {
	 
				if (m.x>x+width-16 && m.x<x+width && m.y>y && m.y<y+16) {
					open ^=1;
					return true;
				}
				if (m.y>y && m.y<y+theight) {
					mouse = Point(m.x-x, m.y-y);
	 
					moving = 1;
					return true;
				}
	 
				return true;
	 
			}
			else if (m.x>x && m.x<x+width && m.y>y+theight && m.y<y+height && open) {
	 
				int mx = m.x-x;
				int my = m.y-y-theight;*/
	std::vector<ToolButton*>::const_iterator it = buttons.begin();
	for (;it!=buttons.end();it++) {

		if ((*it)->contains(mx,my-theight-4)) {
			//			hiliteButton = (*it);
			Message((*it)->id);
			return true;
		}
	}

	//	cout << mx << " " << my << endl;

	std::vector<ToolItem*>::iterator it2 = items.begin();
	for (;it2!=items.end();it2++) {
		if ((*it2)->TakesKey()) {
			if ((*it2)->contains(mx,my)) {
				//			hiliteItem = (*it2);

				(*it2)->MouseClick(button,type,Point<int>(mx,my));
				activeItem = it2;
				hasactiveitem=true;
				return true;
			}
		}
	}


	return true;
	/*
			}
			else {
				moving = 0;
			}
		}
	 
	 
	 
		return false;
	}*/
}

bool ToolBox::KeyboardFunc(SDL_Keysym &keysym) {

	SDL_Keycode key = keysym.sym;

	if (key == SDLK_TAB) {
	
		if (activeItem!=items.end()) {
	//		cout << "TAB1" << endl;
			activeItem++;
			if (activeItem==items.end())
				activeItem = items.begin();

			//cout << "TAB00" << (*activeItem)->rect.text << endl;

			return true;
		} else {
			activeItem == items.begin();
			hasactiveitem=true;
			//cout << "TAB0220" << endl;
			return true;

		}

	} else {
		if (activeItem!=items.end()) {
			(*activeItem)->Key(keysym);
		}
	}

	return true;
}


void ToolBox::addButton(int id_, int image_) {
	int sz = buttons.size();
	buttons.push_back(new ToolButton(id_,image_,2+sz*20,-21));

}

void ToolBox::addStateButton(int &s, int id_, int image_, int modes_) {
	int sz = buttons.size();
	buttons.push_back(new ToolStateButton(s, id_,image_,2+sz*20,-21,modes_));
}
