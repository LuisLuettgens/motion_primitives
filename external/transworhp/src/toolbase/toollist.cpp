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

#include "toollist.h"
#include "../glbase/font.h"



template<class T>
ToolList<T>::ToolList(const ToolRect &rect_, T &r, std::string *Map)
		: ToolData<T>(rect_,r), map(Map), opened(0), openRect(rect_), closeRect(rect_) {

	deleting = 0;

	size = 0;
	index = 0;
	int dd = this->rect.height;
	for (int i=0;map[i][0];i++) {
		openRect.height += dd;
		size++;
	}

}

template<class T>
ToolList<T>::ToolList(const ToolRect &rect_, T &r, const std::vector<std::string> &lst)
		: ToolData<T>(rect_,r), opened(0), openRect(rect_), closeRect(rect_) {

	deleting = 1;
	map = new std::string[lst.size()+1];
	int i=0;
	std::vector<std::string>::const_iterator it = lst.begin();
	for (;it!=lst.end();it++,i++) {
		map[i] = *it;
	}
	map[i] = "";


	size = 0;
	index = 0;
	int dd = this->rect.height;
	for (int i=0;map[i][0];i++) {
		openRect.height += dd;
		size++;
	}

}

template<class T>
ToolList<T>::~ToolList() {
	if (deleting)
		delete []map;
}


template<>
int ToolList<int>::ToValueType(int i) const {
	return i;
}

template<>
color4 ToolList<color4>::ToValueType(int i) const {
	if (map)
	return map[i];
	
	return color4(0);
}

template<>
void ToolList<int>::drawItem(ToolRect &r, int i, double b, float addz) {
	r.text = map[i];
	r.Draw(b,addz);
}

template<>
void ToolList<color4>::drawItem(ToolRect &r, color4 i, double b, float addz) {

	r.text = "";
	r.Draw(b,addz);


	i();
	glBegin(GL_QUADS);

	int dd = 20;
	int p = 3;
	glVertex3f(r.x+p,-r.y-r.height+p,addz+.35);
	glVertex3f(r.x+dd-p,-r.y-r.height+p,addz+.35);
	glVertex3f(r.x+dd-p,-r.y-p,addz+.35);
	glVertex3f(r.x+p,-r.y-p,addz+.35);
	glEnd();


	Tool::colors[5].Dim(.8)();
	Tool::font->printString(ToString(i).c_str(), r.x+r.padding+dd,
	                        -r.y-r.height+r.padding, .3+addz,r.width);


	/*r.text = ToString(i);
	r.x +=dd;
	r.width -=dd;
	//r.Draw(addz);
	r.x -=dd;
	r.width +=dd;
	*/

	/* r.text = std::string("        ") + ToString(i);
	 r.DrawRect();

	 {
	  Rectangle r2 = r;
	  r2.l += 5;
	  r2.t -= 3;
	  r2.b += 3;
	  r2.r = r2.l + 30;
	  r2.z += 0.05f;
	  r2.border=2;
	  r2.SetCol(-1);
	  r2.fill = 1;
	  i();
	  r2.text = "";
	  r2.DrawRect();

	 }
	*/
}

template <class T>
void ToolList<T>::Draw(bool isActive, const Point<int> &mouse) {

	int hilite = this->rect.contains(mouse);

	closeRect.border = Tool::INSET | (hilite?Tool::FILL:0);
	drawItem(closeRect, this->value, .1);


	/*
	 this->back.z = 0.1f;
	 this->back.SetCol(2);
	 this->back.border = 0;
	 this->back.mode = 0;
	 this->back.centertext=0;

	 int anycont=0;
	 if (this->back.Contains(p)) {
	  this->back.mode = 1;
	  anycont = 1;
	 }

	 DrawItem(this->back,this->value);
	*/

	if (!isActive)
		Close();


	if (opened) {
		// this->openRect = this->rect;
		// openRect.SetCol(3);
		// openRect.z +=1.05;

		int dd = closeRect.height;
		/*  openRect.y+= dd;
		  for (int i=0;map[i][0];i++) {
		   openRect.height += dd;
		   openRect.y+= dd;
		  }
		*/
		ToolRect rr = this->openRect;

		rr.border=Tool::SOLID;
		rr.text="";
		rr.y +=dd;
		rr.height -=dd;
		rr.Draw(.4,.4);


		rr = this->closeRect;
		rr.border=Tool::SOLID;		//rr.SetCol(3);
		//rr.z+=1.05;

		for (int i=0;map[i][0];i++) {
			/*rr.mode=0;
			rr.fill=0;
			rr.border = 0;*/
			rr.y+=dd;
			//rr.b-=20;
			//rr.z+=.05;
			// if (rr.Contains(p)) {
			//  rr.mode = 1;
			//rr.fill = 1;
			//  rr.border = 1;
			//  anycont = 1;
			// }
			drawItem(rr,ToValueType(i),0.,.4);
		}
	}
	/* if (anycont==0)
	  opened=0;*/
}



template <class T>
void ToolList<T>::Open() {

	opened=1;
	this->rect = openRect;
}

template <class T>
void ToolList<T>::Close() {

	opened=0;
	this->rect = closeRect;
}

template <class T>
void ToolList<T>::MouseClick(int button, int state,  const Point<int> &p) {

	if (button==1 && state==SDL_MOUSEBUTTONDOWN) {
		if (opened==0)
			Open();

		else {
			ToolRect rr = this->closeRect;
			for (int i=0;map[i][0];i++) {
				rr.y+=closeRect.height;
				if (rr.contains(p.x,p.y)) {
					SetValue(i);

				}
			}
			Close();
		}
	}

	if (button==4 && state==SDL_MOUSEBUTTONDOWN) {
		int a = index-1;
		if (a>=0) {
			SetValue(a);

		}
	}

	if (button==5 && state==SDL_MOUSEBUTTONDOWN) {
		int a = index+1;
		int i=0;
		for (;map[i][0];i++) {}

		if (a<i) {
			SetValue(a);

		}
	}

}

template <class T>
void ToolList<T>::Key(SDL_Keysym &keysym) {

	SDL_Keycode key = keysym.sym;
	// SDLMod mod = keysym.mod;

	//std::cout << (int)c << std::endl;

	switch (key) {

	case 13: // RETURN
		Close();
		break;

	case 27: // ESC
		// editmode = 0;
		break;

	case 127: // ENTF

	case SDLK_UP: {
			int a = index-1;
			if (a<0)
				a=0;
			SetValue(a);

			break;
		}

	case SDLK_DOWN: {
			int a = index+1;
			if (a>size-1)
				a=size-1;
			SetValue(a);

			break;
		}
	default:
		break;

	}
}



template<>
void ToolList<int>::SetValue(int t) {
	value=t;
	index = t;
	Apply();
	if (s) {
		s();
	}
}

template<>
void ToolList<color4>::SetValue(int t) {
	value = color4(map[t]);
	index = t;
	//std::cout << t << ToString(value)<<std::endl;
	Apply();
	if (s) {
		s();
	}
}

template <>
int ToolList<int>::GetValue() {
	return value;
}

template <>
int ToolList<color4>::GetValue() {

	for (int i=0;;i++) {
		if (map[i][0]==0)
			return -1;
		if (ToString(value)==(map[i])) {
			return i;
		}
	}
}


template class ToolList<int>
;
template class ToolList<color4>
;

