#include "toolmenuentry.h"

#include "../base/vectortools.h"
#include "tool.h"

using namespace std;

//void Message(int id);
//void Message(int id, const void *param);

namespace tw {

ToolMenuEntry::ToolMenuEntry(const std::string &t, int i, const std::string &hot, int *sel)
		: text(t), hottext(hot), id(i), pos(0), w2(0), linefrom(0), lineto(0), selected(sel) {

	openentry2 = 0;

	// Automatisch Altkey erkennen
	string::size_type j = text.find('&');
	int alt = ' ';
	if (j!=string::npos) {
		alt = text[j+1];
		string s(1,(char)17);
		text.erase(j,1);

		linefrom = Tool::font->StringLength(string(text,0,j).c_str())-1;
		lineto = Tool::font->StringLength(string(text,0,j+1).c_str())-1;

	}

	w  = Tool::font->StringLength(text.c_str()) + 7;
	w3 = Tool::font->StringLength(hottext.c_str()) ;

	altkey = Tool::toKey(alt);

}

int ToolMenuEntry::GetMenuWidth() {

	int wx = w2+8 + w3 + 13;
	if (w3)
		wx+=20;

	return wx;
}

int ToolMenuEntry::GetMenuHeight() {
	int wy = 0;

	std::vector<ToolMenuEntry*>::const_iterator it = me.begin();
	for (;it!=me.end();it++) {
		if ((*it)->id==-1)
			wy+=4;
		else
			wy +=20;
	}
	return wy;
}


ToolMenuEntry::~ToolMenuEntry() {

	if (me.size())
		ClearPointer(me);

}


void ToolMenuEntry::drawText(int x, int y, float z, int ww, int totalwid, int open,Point<int> &p1, Point<int> &p2) {

	if (open) Tool::drawRect(p1,p2,z-.05,Tool::FILL2 | Tool::INSET);

	Tool::font->printString(text.c_str(),pos+x,-y,z);

	if (hottext!="") {
		Tool::font->printString(hottext.c_str(),x+ww+pos,-y,z);

	}
	if (linefrom<lineto) {
		glBegin(GL_LINES);
		glVertex3f(pos+x+linefrom, -y-1, z+.1);
		glVertex3f(pos+x+lineto, -y-1, z+.1);
		glEnd();
	}

	if (me.size() && totalwid) {
		char buf[] = {static_cast<char>(FFont::GORIGHT),0};
		Tool::font->printString(buf,x+pos+totalwid-6,-y,z);

	}

	if (selected && totalwid) {
		char buf[] = {static_cast<char>(FFont::UNCHECK),0};
		if (*selected)
			buf[0] = FFont::CHECK;
		Tool::font->printString(buf,pos+x-14,-y,z);
	}


}


void ToolMenuEntry::DrawOpenMenu(const ToolMenu *m, int x, int y, float z) {

	int wx = GetMenuWidth();
	int wy = GetMenuHeight();

	Point<int> p1(x,    y);
	Point<int> p2(wx+x, y+wy);

	Tool::drawRect(p1,p2,z,Tool::FILL | Tool::OUTSET);

	std::vector<ToolMenuEntry*>::const_iterator it = me.begin();
	for (int i=18;it!=me.end();it++,i+=20) {

		if ((*it)->id==-1) {
			Point<int> p1(x+4,    y+i-16);
			Point<int> p2(x+wx-4, y+i-15);
			Tool::drawRect(p1,p2,z+.1, Tool::FILL2 | Tool::INSET);
			i-=16;
		}
		else {

// 			if (openentry2 == *it) {
 				Point<int> p1(x+2,    y+i-16);
				Point<int> p2(x+wx-2, y+i);
//
// 				Tool::drawRect(p1,p2,z+.05,Tool::FILL2 | Tool::INSET);
// 			}
			glColor4fv(Tool::colors[4].GetData());
			(*it)->drawText(x+ 5+13,y+i-4,z+.1,w2+20,w2+w3+20, openentry2 == *it,p1, p2);

			if (openentry2 == *it) {
				if (openentry2->me.size()) {
					yy = y+i-18+4;
					openentry2->DrawOpenMenu(m,x+wx-4,yy,z+.2);
				}
			}
		}
	}
}


bool ToolMenuEntry::MouseMotion(const Point<int> &mouse, int x, int y) {

	int wx = GetMenuWidth();
	int wy = GetMenuHeight();

	Point<int> p1(x,    y);
	Point<int> p2(x+wx, y+wy);

	if (mouse.isin(p1,p2)) {

		std::vector<ToolMenuEntry*>::iterator it = me.begin();
		for (int i=18;it!=me.end();it++,i+=20) {

			if ((*it)->id==-1) {
				i-=16;
			}
			else {
				Point<int> p1(x,    y+i-20);
				Point<int> p2(x+wx, y+i);

				if (mouse.isin(p1,p2)) {
					openentry2 = *it;
					openentry2->openentry2 = 0;
					break;
				}
			}
		}

		return true;
	}

	if (openentry2) {
		if (openentry2->me.size()) {
			if (openentry2->MouseMotion(mouse,x+wx-4,yy))
				return true;
		}
	}



	openentry2=0;

	return false;
}


bool ToolMenuEntry::MouseButton(const Point<int> &mouse, int x, int y, Uint8 type, Uint8 button) {

	int wx = GetMenuWidth();
	int wy = GetMenuHeight();

	Point<int> p1(x,    y);
	Point<int> p2(x+wx, y+wy);

	if (openentry2) {
		if (openentry2->me.size()) {
			if (openentry2->MouseButton(mouse,x+wx-4,y,type,button))
				return true;
		}
	}

// MATTHIAS!!!
/*	if (openentry2) {
		if (!openentry2->openentry2) {
			if (button == SDL_BUTTON_WHEELUP && type == SDL_MOUSEBUTTONUP) {
				KeyBoardFunc(SDLK_UP);
				return true;
			}

			if (button == SDL_BUTTON_WHEELDOWN && type == SDL_MOUSEBUTTONUP) {
				KeyBoardFunc(SDLK_DOWN);
				return true;
			}
		}
	}*/



	if (mouse.isin(p1,p2)) {

		std::vector<ToolMenuEntry*>::iterator it = me.begin();
		for (int i=18;it!=me.end();it++,i+=20) {

			if ((*it)->id==-1) {
				i-=16;
			}
			else {
				p1.y = y+i-20;
				p2.y = y+i;

				if ( mouse.isin(p1,p2) && type == 1 /* up */ && button == SDL_BUTTON_LEFT ) {
					if ((*it)->MouseButton1(mouse, x, y, type, button)) {
						openentry2 = 0;
						return true;
					}
				}
			}
		}

		//return true;
	}

	return false;
}


void ToolMenuEntry::ResetSelection() {
	if (me.size()) {
		openentry2 = *me.begin();
		openentry2->openentry2 = 0;
	}
	else
		openentry2 = 0;

}



bool ToolMenuEntry::KeyBoardFunc(SDL_Keycode key) {

	if (openentry2)
		if (openentry2->KeyBoardFunc(key))
			return true;


	if (key == SDLK_LEFT) {
		if (me.size() && openentry2) {
			openentry2 = 0;
			return true;
		}


	}

	if (key == SDLK_RIGHT) {
		if (me.size() && openentry2==0) {
			openentry2 = *me.begin();
			openentry2->openentry2=0;
			return true;
		}

	}

	if (key == SDLK_DOWN) {
		if (me.size() && openentry2) {

			std::vector<ToolMenuEntry*>::iterator it = me.begin();
			for (;it!=me.end();it++) {

				if (openentry2==*it) {
					do {
						it++;
						if (it==me.end())
							it = me.begin();

						openentry2 = *it;
						openentry2->openentry2=0;

					}
					while (openentry2->id == -1);

					break;
				}
			}
			return true;
		}
	}

	if (key == SDLK_UP) {
		if (me.size() && openentry2) {

			std::vector<ToolMenuEntry*>::reverse_iterator it = me.rbegin();
			for (;it!=me.rend();it++) {

				if (openentry2==*it) {
					do {
						it++;
						if (it==me.rend())
							it = me.rbegin();

						openentry2 = *it;
						openentry2->openentry2=0;

					}
					while (openentry2->id == -1);

					break;
				}
			}
			return true;
		}
	}

	if (key == SDLK_RETURN) {

		if (id) {
			Message(id);
			return true;
		}

	}


	return false;
}











ToolMenuSliderEntry::ToolMenuSliderEntry(const std::string &t, int i, const std::string &hot, float *val_, float minval_, float maxval_)
		: ToolMenuEntry(t,i,hot,0),
		val(val_), minval(minval_), maxval(maxval_) {

	openentry2 = 0;

	// Automatisch Altkey erkennen
	string::size_type j = text.find('&');
	int alt = ' ';
	if (j!=string::npos) {
		alt = text[j+1];
		string s(1,(char)17);
		text.erase(j,1);

		linefrom = Tool::font->StringLength(string(text,0,j).c_str())-1;
		lineto = Tool::font->StringLength(string(text,0,j+1).c_str())-1;

	}

	w  = Tool::font->StringLength(text.c_str()) + 7;
	w3 = Tool::font->StringLength(hottext.c_str());

	w += 100;
	w2 = 100;


	altkey = Tool::toKey(alt);

}


bool ToolMenuSliderEntry::KeyBoardFunc(SDL_Keycode key) {

	if (key == SDLK_PLUS) {
		(*val)++;
		if (*val>maxval) *val=maxval;
	}

	if (key == SDLK_MINUS) {
		(*val)--;
		if (*val<minval) *val=minval;
	}


	return ToolMenuEntry::KeyBoardFunc( key);

}


void ToolMenuSliderEntry::drawText(int x, int y, float z, int ww, int totalwid, int open, Point<int> &p1, Point<int> &p2) {

	Point<int> p3(x+w2-4,p2.y);

	if (open) Tool::drawRect(p1,p3,z-.05,Tool::FILL2 | Tool::INSET);



	Tool::font->printString(text.c_str(),pos+x,-y,z);

	if (hottext!="") {
		Tool::font->printString(hottext.c_str(),x+ww+pos,-y,z);

	}
	if (linefrom<lineto) {
		glBegin(GL_LINES);
		glVertex3f(pos+x+linefrom, -y-1, z+.1);
		glVertex3f(pos+x+lineto, -y-1, z+.1);
		glEnd();
	}

	if (me.size() && totalwid) {
		char buf[] = {static_cast<char>(FFont::GORIGHT),0};
		Tool::font->printString(buf,x+pos+totalwid-6,-y,z);
	}

	if (totalwid) {
		//int wx = GetMenuWidth();

		Point<int> p1a(x+w2, p1.y);
		Point<int> p2a(x+w2+72, p2.y);

		Tool::drawRect(p1a,p2a,z+.005, Tool::FILL2 | Tool::OUTSET);

		char buf[20];
		sprintf(buf,"%10.2f",*val);
		Tool::font->printString(buf,pos+x+w2+10,-y,z+.01);
	}

}



bool ToolMenuEntry::MouseButton1(const Point<int> &/*mouse*/, int /*x*/, int /*y*/, Uint8 /*type*/, Uint8 /*button*/) {

	if (id)
		Message(id);

	return true;
}


bool ToolMenuSliderEntry::MouseButton1(const Point<int> &/*mouse*/, int /*x*/, int /*y*/, Uint8 /*type*/, Uint8 /*button*/) {

	if (id)
		Message(id);

	cout << "IN" << endl;

	return false;
}

}
