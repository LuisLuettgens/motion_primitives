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
#include "toolitem.h"

#include "GL/gl.h"
#include "conversion.h"

ToolItem::ToolItem() : rect(0,0,0,0), s(0) {}

ToolItem::ToolItem(const ToolRect &rect_)
		: rect(rect_), s(0) {}

ToolItem::~ToolItem() {}

bool ToolItem::contains(int x_,int y_) {
	return rect.contains(x_,y_);
}

void ToolItem::Draw(bool, const Point<int>&) {}

void ToolItem::MouseClick(int, int, const Point<int>&) {}
void ToolItem::Key(SDL_Keysym&) {}
//void SpecialKey(Uint8*) {}
bool ToolItem::TakesKey() {
	return false;
}
void ToolItem::Apply() {}
void ToolItem::Discard() {}

void ToolItem::Signal(signal1 s_) {
	s=s_;
}


template <class T>
void ToolLabel<T>::Draw(bool isActive, const Point<int> &mouse) {

	Update();

	rect.Draw(0);

}


template <class T>
ToolLabel<T>::ToolLabel(const ToolRect &rect_, T &r)
		: ToolItem(rect_), ref(r) {}

template <>
ToolLabel<std::string>::ToolLabel(const ToolRect &rect_, std::string &r)
		: ToolItem(rect_), ref(r) {}

template <>
ToolLabel<const std::string>::ToolLabel(const ToolRect &rect_, const std::string &r)
		: ToolItem(rect_), ref(r) {
	rect.text=r;
}

template <class T>
void ToolLabel<T>::Update() {
	rect.text = ToString(this->ref);
}

template<>
void ToolLabel<std::string>::Update() {
	rect.text = this->ref;
}

template<>
void ToolLabel<const std::string>::Update() {}


template class ToolLabel<double>
;
template class ToolLabel<std::string>
;
template class ToolLabel<int>
;
template class ToolLabel<const std::string>
;









void drawRect2(const Point<int> &p1,const Point<int> &p2,float z, int mode, int a, int b) {

	if (mode & (1+2+Tool::FILLDARK+Tool::FILLHI)) {
		if (mode&2)
			Tool::colors[0]();
		else if (mode&Tool::FILLDARK)
			Tool::colors[7]();
		else if (mode&Tool::FILLHI)
			Tool::colors[7]();
		else
			Tool::colors[2]();

		glBegin(GL_QUADS);

		glVertex3f(p1.x, -p1.y, z);
		glVertex3f(p1.x, -p2.y, z);

		glVertex3f(p2.x, -p2.y, z);
		glVertex3f(p2.x, -p1.y, z);
		glEnd();
	}

	if (mode & (4+8+16)) {
		Tool::colors[mode&8?3:1]();
		glBegin(GL_LINES);
		glVertex3f(p2.x+1, -p1.y, z+.1);
		glVertex3f(p1.x+b, -p1.y, z+.1);

		glVertex3f(p1.x+a, -p1.y, z+.1);
		glVertex3f(p1.x, -p1.y, z+.1);

		glVertex3f(p1.x, -p1.y, z+.1);
		glVertex3f(p1.x, -p2.y, z+.1);

		Tool::colors[mode&4?3:1]();
		glVertex3f(p1.x, -p2.y, z+.1);
		glVertex3f(p2.x, -p2.y, z+.1);

		glVertex3f(p2.x, -p2.y, z+.1);
		glVertex3f(p2.x, -p1.y, z+.1);
		glEnd();
	}
}

template <class T>
void ToolFrame<T>::Draw(bool isActive, const Point<int> &mouse) {

	Update();


	int ww = Tool::font->StringLength(r2.text.c_str());

	drawRect2(Point<int>(rect.x,rect.y),Point<int>(rect.x+rect.width,rect.y+rect.height),
	          0.3,Tool::OUTSET,6,8+ww);

	
	drawRect2(Point<int>(rect.x-1,rect.y-1),Point<int>(rect.x+rect.width+1,rect.y+rect.height+1),
	          0.3,Tool::INSET,6,8+ww);


	//r2.border=Tool::FILLD;
	r2.Draw(0,.9);
}


template <class T>
ToolFrame<T>::ToolFrame(const ToolRect &rect_, T &r)
		: ToolItem(rect_), ref(r), r2(rect_) {

	rect.y+=12;
	r2.y+=12;
	rect.height-=12;
	r2.height-=12;

	r2.x+=4;
	r2.y=rect_.y+rect.height;
	r2.height=16;

}

template <>
ToolFrame<std::string>::ToolFrame(const ToolRect &rect_, std::string &r)
		: ToolItem(rect_), ref(r), r2(rect_) {

	rect.y+=12;
	r2.y+=12;
	rect.height-=12;
	r2.height-=12;

	r2.y=rect_.y+rect.height;
	r2.x+=4;
	r2.height=16;

}

template <>
ToolFrame<const std::string>::ToolFrame(const ToolRect &rect_, const std::string &r)
		: ToolItem(rect_), ref(r), r2(rect_) {

	rect.y+=12;
	r2.y+=12;
	rect.height-=12;
	r2.height-=12;

	r2.text = r;
	r2.y=rect_.y;//+rect.height;
	r2.x+=4;
	r2.height=16;

}

template <class T>
void ToolFrame<T>::Update() {
	r2.text = ToString(this->ref);
}

template<>
void ToolFrame<std::string>::Update() {
	r2.text = this->ref;
}

template<>
void ToolFrame<const std::string>::Update() {}




template class ToolFrame<std::string>
;

template class ToolFrame<const std::string>
;


