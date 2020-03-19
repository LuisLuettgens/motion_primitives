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

#define _USE_MATH_DEFINES
#include <cmath>

#include "toollistview.h"

#include "GL/gl.h"
#include "conversion.h"
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include "../base/twstatus.h"


using namespace std;


ToolListview::ToolListview(const std::string &t, int w, int h, bool write) : ToolWindow(t) {

	message_id = 0;
	moving = 0;
	width = 20;
	height = 140;
	x = (w - width)/2;
	y = (h - height)/2;

	open = 1;
	theight = 18;

	AddIcon(ToolIcon::CLOSE);
	if (write)
		AddIcon(ToolIcon::WRITETXT);

}

void ToolListview::AddColumn(const std::string &s, int w) {
	titles.push_back(new ToolListcolumn(s,w));
	width += w;
	x -= w/2;
}


ToolListview::~ToolListview() {
	// ClearPointer(buttons);
	// ClearPointer(items);
	clearpointer(list);
	clearpointer(titles);
}






bool ToolListview::mousebutton(int mx, int my, int button, int type) {

	if (button==1) {
		if (mx>=10 && mx <= width-10) {
			int index = (my - 10) / 16;

			if (message_id) {
				if (index>=1) {
					ToolListdata *tld = list[index-1];
					ToolListitem *tli = tld->data[0];
					Message(message_id,tli->GetString().c_str());
				}
			} else {
				bool lastopen = false;

				std::vector<ToolListdata*>::iterator it3 = list.begin();
				for (;it3!=list.end();it3++) {

					if ((*it3)->level == 0) {
						lastopen = (*it3)->open;
					}

					if ((*it3)->level ==0 || lastopen) {
						index--;
					}

					if (index==0) {
						// cout << (*(*it3)->data.begin())->GetString() << endl;
						(*it3)->open = !(*it3)->open;
						break;
					}
				}

			}


		}
	}
	return true;

}

//void ToolListview::WriteImage() {gl->iw->WriteTool(this);}

void ToolListview::WriteText() {

	stringstream a;
	a << Tool::filename_prefix + string("_") + title + ".dat";


	ofstream of(a.str().c_str());

	stringstream st;
	st << a.str();
	MyStatus("Writing", st.str(), Status::NORMAL);

	of << left;

	{
		vector<ToolListcolumn*>::iterator it = titles.begin();
		for (;it!=titles.end();it++) {

			of << setw(20) << (*it)->text;
		}
	}
	of << endl;


	std::vector<ToolListdata*>::iterator it3 = list.begin();
	for (;it3!=list.end();it3++) {

		std::vector<ToolListitem*>::iterator it2 = (*it3)->data.begin();
		for (;it2!=(*it3)->data.end();it2++) {
			of << setw(20) << (*it2)->GetString();
		}

		of << endl;
	}


}


void ToolListview::Draw(bool isActive, int HT, const Point<int> &mouse, bool border) {

	ToolWindow::Draw(isActive, HT, mouse, border);

	innerProj(HT);

	int xx = 10, yy = 10;

	ToolRect r(xx,yy,0, 15);
	std::vector<ToolListcolumn*>::iterator it = titles.begin();
	for (;it!=titles.end();it++) {

		r.width = (*it)->width-1;
		r.text = (*it)->text;
		r.border = Tool::OUTSET;
		r.Draw(.3);

		r.x += (*it)->width;
	}


	r = ToolRect(xx,yy,0, 15);


	bool lastopen = false;

	std::vector<ToolListdata*>::iterator it3 = list.begin();
	for (;it3!=list.end();it3++) {

		if ((*it3)->level == 0) {
			lastopen = (*it3)->open;
		}

		if ((*it3)->level == 0 || lastopen) {
			r.y += 16;
			r.x = xx;

			it = titles.begin();
			std::vector<ToolListitem*>::iterator it2 = (*it3)->data.begin();
			for (;it2!=(*it3)->data.end();it2++,it++) {
				r.width = (*it)->width-1;
				r.border = Tool::INSET;

				// r.text = (*it2)->ToString();
				(*it2)->Draw(r);
				r.x += (*it)->width;
			}
		}
	}


	scrollmax = r.y + 40;
	if (scrollmax>height) {
		if (scrollcur > scrollmax - height)
			scrollcur = scrollmax-height;
	}
	/*
	// Tiefe von 0 bis 1

	// Titelzeile
	 glBegin(GL_QUADS);

	 toolcolors[0].dim(moving?.1:.0);

	 glVertex3f(x,-y,.0);
	 glVertex3f(x,-y-theight,.0);

	 toolcolors[0].dim(moving?.15:.05);

	 glVertex3f(x+width,-y-theight,.0);
	 glVertex3f(x+width,-y,.0);


	// Textfeld

	 if (open) {
	  toolcolors[1].dim(.0);

	  glVertex3f(x,-y-theight,.0);
	  glVertex3f(x,-y-height,.0);

	  toolcolors[1].dim(.05);

	  glVertex3f(x+width,-y-height,.0);
	  glVertex3f(x+width,-y-theight,.0);

	 }
	 glEnd();

	 if (isActive)
	  glLineWidth(2);
	 else {
	  glLineWidth(1);
	 }

	 glDisable(GL_LINE_SMOOTH);

	 toolcolors[1].dim(.7);
	 glBegin(GL_LINE_STRIP);
	 glVertex3f(x,-y,.1);
	 glVertex3f(x+width+.1,-y,.1);
	 if (open) {
	  glVertex3f(x+width+.1,-y-height,.1);
	  glVertex3f(x,-y-height,.1);
	 }
	 else {
	  glVertex3f(x+width+.1,-y-theight,.1);
	  glVertex3f(x,-y-theight,.1);
	 }
	 glVertex3f(x,-y,.1);
	 glEnd();

	*/
	/*

	 glTranslatef(x,-y-theight,0);

	 toolcolors[2].alpha(.8);

	 int w = gl->font->StringLength(title.c_str());
	 gl->font->printString(title.c_str(), (width-w)/2, 3, .1);


	// Icon rechts
	 int icon = 16;
	 float tx1 = open?0:1./8.;
	 float ty1 = 0.;
	 float tx2 = open?1./8.:2./8.;
	 float ty2 = 1./8.;
	 glEnable(GL_TEXTURE_2D);
	 gl->th.Bind(7);
	 glBegin(GL_QUADS);

	 toolcolors[0].dim(moving?.3:.2);

	 glTexCoord2f(tx1,ty1);
	 glVertex3f(width-icon,icon,.5);
	 glTexCoord2f(tx1,ty2);
	 glVertex3f(width-icon,0,.5);
	 glTexCoord2f(tx2,ty2);
	 glVertex3f(width,0,.5);
	 glTexCoord2f(tx2,ty1);
	 glVertex3f(width,icon,.5);
	 glEnd();
	 glDisable(GL_TEXTURE_2D);


	// gl->th.Bind(7);

	// Iconleiste

	 glLineWidth(1);


	 if (open) {
	  std::vector<ToolButton*>::const_iterator it = buttons.begin();
	  int i=0;
	  for (;it!=buttons.end();it++,i++) {
	   toolcolors[0].dim(moving?.3:.2);
	   (*it)->Draw();
	  }

	  Update();


	  std::vector<ToolItem*>::const_iterator it2 = items.begin();
	  for (;it2!=items.end();it2++) {
	   (*it2)->Draw(activeItem == *it2 && isActive);
	  }


	  toolcolors[2].alpha(.6);
	  DrawItem();

	 }
	 glLineWidth(1);
	*/
}


/*void ToolListview::Apply() {

 std::vector<ToolItem*>::const_iterator it2 = items.begin();
 for (;it2!=items.end();it2++) {
  (*it2)->Apply();
 }

}*/

/*
bool ToolListview::KeyboardFunc(SDL_keysym &keysym) {

 if (activeItem) {
  activeItem->Key(keysym);
 }

 return true;
}
*/




std::string ToStringPerc(double v) {

	char buf[40];
	sprintf(buf,"%.2f%%",v*100.);
	return std::string(buf);

}

std::string ToStringWinkel(double v) {

	char buf[40];
	sprintf(buf,"%.2f�",v*180./M_PI);
	return std::string(buf);

}


void ToolListdata::AddItem(ToolListitem* ti) {

	data.push_back(ti);

}

template <>
const std::string ToolListD<int>::GetString() const {
	char buf[40];
	sprintf(buf,"%d",this->value);
	return std::string(buf); //String((T) this->value);
}
template <>
const std::string ToolListD<double>::GetString() const {
	if (mode==0) {
		char buf[40];
		sprintf(buf,"%.3f",this->value);
		return std::string(buf); //String((T) this->value);
	}
	if (mode==1)
		return ToStringPerc(this->value);
	if (mode==2)
		return ToStringWinkel(this->value);
	return "";
}
template<>
const std::string ToolListD<std::string>::GetString() const {
	return value;
}
template<>
const std::string ToolListD<const std::string>::GetString() const {
	return (value);
}

template<class T>
void ToolListD<T>::Draw(ToolRect &r) {
	r.text = GetString();
	r.Draw();

}



template <>
const std::string ToolListR<int>::GetString() const {
	char buf[40];
	sprintf(buf,"%d",this->ref);
	return std::string(buf); //String((T) this->value);
}
template <>
const std::string ToolListR<double>::GetString() const {
	if (mode==0) {
		char buf[40];
		sprintf(buf,"%.3f",this->ref);
		return std::string(buf); //String((T) this->value);
	}
	if (mode==1)
		return ToStringPerc(this->ref);
	if (mode==2)
		return ToStringWinkel(this->ref);
	return "";
}

template<>
const std::string ToolListR<std::string>::GetString() const {
	return ref;
}
template<>
const std::string ToolListR<const std::string>::GetString() const {
	return (ref);
}

template<class T>
void ToolListR<T>::Draw(ToolRect &r) {
	r.text = GetString();
	r.Draw();

}




template<class T>
void ToolListC<T>::Draw(ToolRect &r) {

	col.Dim(0)();
	glBegin(GL_QUADS);

	int dd = 20;
	int p = 3;
	glVertex3f(r.x+p,-r.y-r.height+p,.25);
	glVertex3f(r.x+dd-p,-r.y-r.height+p,.25);
	glVertex3f(r.x+dd-p,-r.y-p,.25);
	glVertex3f(r.x+p,-r.y-p,.25);
	glEnd();

	r.text = "";
//	r.x +=dd;
//	r.width -=dd;
	r.Draw();
//	r.x -=dd;

	ToolRect r2(r);
	r2.text = this->GetString();
	r2.x +=dd;
	r2.width -=dd;
	r2.border = 0;
	r2.Draw();


}

template class ToolListD<int>
;
template class ToolListD<double>
;
template class ToolListD<const std::string>
;

template class ToolListR<int>
;
template class ToolListR<double>
;
template class ToolListR<const std::string>
;

template class ToolListC<const std::string>
;
