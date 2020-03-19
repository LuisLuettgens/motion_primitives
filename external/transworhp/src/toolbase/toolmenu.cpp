#include "toolmenu.h"

#include "../base/vectortools.h"
#include "toolmenuentry.h"

#include "toolrect.h"

using namespace std;

//void Message(int id);
//void Message(int id, const void *param);

namespace tw {

ToolMenu::ToolMenu()
		: openentry(0), opened(0), visible(1),yy(100), menuposct(5) {}


ToolMenuEntry *ToolMenu::AddMenu(const std::string &s, int id) {

	ToolMenuEntry *m = new ToolMenuEntry(s,id,"",0);

	m->pos = menuposct;
	menuposct += m->w+1;

	me.push_back(m);
	return m;
}

ToolMenuEntry *ToolMenu::AddSeparator(ToolMenuEntry *th) {

	ToolMenuEntry *m = new ToolMenuEntry("",-1,"",0);

	/*if (m->w > th->w2)
		th->w2 = m->w;*/

	th->me.push_back(m);

	return m;
}
/*
MenuEntry *Menu::AddMenu(MenuEntry *th, const std::string &s, int id) {

	MenuEntry *m = new MenuEntry(s,id,"",0);

	if (m->w > th->w2)
		th->w2 = m->w;

	th->me.push_back(m);

	return m;
}*/

ToolMenuEntry *ToolMenu::AddMenu(ToolMenuEntry *th, const std::string &s, int id, const std::string &hot, int *selected) {

	if (hot!="")
		mk.push_back(ToolMenuKey(hot,id));

	ToolMenuEntry *m = new ToolMenuEntry(s,id,hot,selected);

	if (m->w > th->w2)
		th->w2 = m->w;

	if (hot!="") {
		if (m->w3 > th->w3)
			th->w3 = m->w3;
	}
	th->me.push_back(m);

	return m;
}

ToolMenuEntry *ToolMenu::AddMenu(ToolMenuEntry *th, const std::string &s, int id, float *val, float minval, float maxval,const std::string &hot) {

	if (hot!="")
		mk.push_back(ToolMenuKey(hot,id));

	ToolMenuEntry *m = new ToolMenuSliderEntry(s,id,hot,val,minval,maxval);

	if (m->w > th->w2)
		th->w2 = m->w;

	if (hot!="") {
		if (m->w3 > th->w3)
			th->w3 = m->w3;
	}
	th->me.push_back(m);

	return m;
}


ToolMenuEntry *ToolMenu::AddMenu(ToolMenuEntry *th, const std::string &s, int id, SmoothMovement &sm,const std::string &hot) {

	if (hot!="")
		mk.push_back(ToolMenuKey(hot,id));

	ToolMenuEntry *m = new ToolMenuSliderEntry(s,id,hot,&sm.value,sm.minval,sm.maxval);

	if (m->w > th->w2)
		th->w2 = m->w;

	if (hot!="") {
		if (m->w3 > th->w3)
			th->w3 = m->w3;
	}
	th->me.push_back(m);

	return m;
}

ToolMenu::~ToolMenu() {

	ClearPointer(me);

}


void ToolMenu::Draw(int w) const {

	if (!visible) {
		return;
	}

	float alpha = 1.0f;
	if (visible < 20) {
		alpha = visible/20.;
	}

	//if (!ToolMenuEntry::font)
	//	return;
	if (!me.size()) {
		return;
	}

	const float z = 11.0f;

	glDisable(GL_LINE_SMOOTH);

	const int w1 = me.back()->pos + me.back()->w;

	glColor4fv(Tool::colors[0].Alpha(alpha).GetData());

	glBegin(GL_QUADS);

	glVertex3f(0, 0, z +.1);
	glVertex3f(0, -20, z+.1);

	glVertex3f(w1, -20, z+.1);
	glVertex3f(w1, 0, z+.1);


	glVertex3f(w1, 0, z+.1);
	glVertex3f(w1, -20, z+.1);

	glColor4fv(Tool::colors[0].Alpha(.0).GetData());
	glVertex3f(w, -20, z+.1);
	glVertex3f(w, 0, z+.1);

	glEnd();

	for (auto it = me.cbegin(); it != me.end(); ++it) {

		if (openentry==*it) {
			Point<int> p1(openentry->pos,2);
			Point<int> p2(openentry->pos+openentry->w,18);
			Tool::drawRect(p1,p2,z+.15,Tool::FILL | Tool::INSET);
		}
		Point<int> p1,p2;
		glColor4fv(Tool::colors[4].Alpha(alpha).GetData());
		(*it)->drawText(5,15,z+.2,0,0,0,p1,p2);

	}

	{
		ToolMenuEntry entry(infotext,0,"",0);
		entry.pos = w-entry.w-5;
		/*{
			Point<int> p1(entry.pos,2);
			Point<int> p2(entry.pos+entry.w,18);
			Tool::drawRect(p1,p2,z+.15,Tool::FILL | Tool::INSET);
		}*/
		Point<int> p1,p2;
		glColor4fv(Tool::colors[5].Alpha(alpha).GetData());
		entry.drawText(5,15,z+.2,0,0,0,p1,p2);
	}

	if (openentry && opened) {
		openentry->DrawOpenMenu(this,openentry->pos,18,z+.2);
	}
}


void ToolMenu::Timer(int t) {

	if (yy>=20) {
		if (!openentry && !opened && visible)
			visible--;
	}

	status.Timer(t);

}


bool ToolMenu::MouseMotion(SDL_MouseMotionEvent &m) {

	if (m.y<20) {
		visible = 100;
	}
	yy = m.y;

	if (!visible)
		return false;

	Point<int> mouse(m.x,m.y);

	std::vector<ToolMenuEntry*>::iterator it = me.begin();
	for (;it!=me.end();it++) {

		Point<int> p1((*it)->pos,2);
		Point<int> p2((*it)->pos+(*it)->w,18);
		if (mouse.isin(p1,p2)) {

			if (openentry!=*it) {
				openentry=*it;
				openentry->openentry2 = 0;
			}
			return true;
		}
	}

	if (openentry && opened) {
		if (openentry->MouseMotion(mouse,openentry->pos,18))
			return true;
	}

	opened = 0;
	openentry = 0;

	return 0;

}


bool ToolMenu::MouseButton(SDL_MouseButtonEvent &m) {

	if (!visible)
		return false;

	Point<int> mouse(m.x,m.y);

	if (openentry && opened) {

		if (openentry->MouseButton(mouse,openentry->pos,18,m.type,m.button)) {

			if (m.type == SDL_MOUSEBUTTONUP && m.button == SDL_BUTTON_LEFT) {
				opened = 0;
				openentry = 0;
			}
			return true;
		}
	}

	if (m.y<20) {
		std::vector<ToolMenuEntry*>::iterator it = me.begin();
		for (;it!=me.end();it++) {

			Point<int> p1((*it)->pos,2);
			Point<int> p2((*it)->pos+(*it)->w,18);
			if (mouse.isin(p1,p2)
			        && m.type == SDL_MOUSEBUTTONDOWN
			        && m.button == SDL_BUTTON_LEFT) {

				opened ^= 1;

				/*if (openentry!=*it)
					openentry= *it;
				else
					openentry=0;
				*/
				return true;
			}
		}
		return true;
	}

	return false;
}




bool ToolMenu::KeyboardFunc(SDL_Keysym &keysym) {

	SDL_Keycode key = keysym.sym;
	Uint16 mod = keysym.mod;

	// Hotkeys checken
	std::vector<ToolMenuKey>::iterator it = mk.begin();
	for (;it!=mk.end();it++) {

		if (it->Check(key,mod))
			return true;


	}


	if (key == SDLK_LALT  || key == SDLK_RALT ) {

		if (openentry) {
			openentry = 0;
		} else {
			openentry = *me.begin();
			opened = 0;
			visible = 100;
		}

		return true;
	}


	if (!openentry)
		return false;


	// Altkeys
	if (!opened) {
		std::vector<ToolMenuEntry*>::iterator it = me.begin();
		for (;it!=me.end();it++) {

			if ((*it)->altkey == key) {
				openentry = *it;
				opened = 1;
				if (openentry->me.size())
					openentry->openentry2 = *openentry->me.begin();
				else
					openentry->openentry2 = 0;

				return true;
			}
		}

	} else {
		std::vector<ToolMenuEntry*>::iterator it = openentry->me.begin();
		for (;it!=openentry->me.end();it++) {

			if ((*it)->altkey == key) {

				openentry->openentry2 = *it;

				if ((*it)->id) {

					Message((*it)->id);
				}

				opened = 0;
				openentry = 0;


				return true;
			}

		}
	}




	if (key == SDLK_LEFT) {

		if (openentry->openentry2)
			if (openentry->openentry2->KeyBoardFunc(key))
				return true;

		std::vector<ToolMenuEntry*>::reverse_iterator it = me.rbegin();
		for (;it!=me.rend();it++) {

			if (openentry==*it) {
				it++;
				if (it==me.rend())
					it = me.rbegin();

				openentry = *it;
				openentry->ResetSelection();
				break;
			}
		}
	}

	if (key == SDLK_RIGHT) {

		if (openentry->openentry2)
			if (openentry->openentry2->KeyBoardFunc(key))
				return true;

		std::vector<ToolMenuEntry*>::iterator it = me.begin();
		for (;it!=me.end();it++) {

			if (openentry==*it) {
				it++;
				if (it==me.end())
					it = me.begin();

				openentry = *it;
				openentry->ResetSelection();
				break;
			}
		}
	}


	if (key == SDLK_DOWN) {

		if (openentry->openentry2)
			if (openentry->openentry2->KeyBoardFunc(key))
				return true;

		if (!opened) {
			opened = 1;
			openentry->ResetSelection();

		} else {

			if (openentry->openentry2==0) {
				openentry->openentry2 = *(openentry->me.begin());
				openentry->openentry2->openentry2=0;
			} else {
				openentry->KeyBoardFunc(key);
			}
		}
	}

	if (key == SDLK_UP) {

		if (openentry->openentry2)
			if (openentry->openentry2->KeyBoardFunc(key))
				return true;

		if (!opened) {
			opened = 1;
			openentry->ResetSelection();
		} else {

			if (openentry->openentry2==0) {
				openentry->openentry2 = *(openentry->me.rbegin());
				openentry->openentry2->openentry2=0;
			} else {
				openentry->KeyBoardFunc(key);

			}
		}
	}

	if (key == SDLK_ESCAPE) {
		opened = 0;
		openentry = 0;
	}

	if (key == SDLK_RETURN) {

		if (openentry->KeyBoardFunc(key)) {
			opened = 0;
			openentry = 0;
		}
	}


	if (opened) {
		if (openentry->openentry2->KeyBoardFunc(key))
				return true;
	}

	return true;

}


bool ToolMenu::SpecialKeys(const Uint8 */*keys*/) {

	if (openentry) {
		return true;
	}

	return false;
}

}
