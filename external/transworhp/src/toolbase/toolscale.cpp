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
#include <string>

#include "GL/gl.h"

#include "conversion.h"
#include "toolscale.h"



template <class T>
std::string ToolScale<T>::GetString() {
	return ToString(this->value);
}

template<>
void ToolScale<int>::SetString(const std::string &s) {
	value = ToInt(s);
}


template <class T>
void ToolScale<T>::MouseClick(int button, int state, const Point<int> &p) {


	if (button==1 && state==SDL_MOUSEBUTTONDOWN) {

		this->value++;
		if (this->value==3)
			this->value = 0;

		this->ref = this->value;

	}





}




template <class T>
void ToolScale<T>::Key(SDL_Keysym &keysym) {

	SDL_Keycode key = keysym.sym;
	// SDLMod mod = keysym.mod;

	//std::cout << (int)c << std::endl;

	switch (key) {
	case 8: // <--

		break;

	case 13: // RETURN
		// editmode = 0;
		break;

	case 27: // ESC
		// editmode = 0;
		break;

	case 32: // leer
		// editmode = 0;
		break;

	case 127: // ENTF

		break;

	default:
		break;
	}
}



template <class T>
int ToolScale<T>::getht(double i) {

//int xx = this->rect.x+this->rect.width-icon;
// int yy = -(this->rect.y);
// int h = this->rect.height-1-10;

//(step-i)/(double)step

	return (int)(-this->rect.y-(this->rect.height-12)*i+this->rect.padding-8);

//yy-h*i/5.+5;

}

template <class T>
void ToolScale<T>::Draw(bool isActive, const Point<int> &mouse) {

	int icon = 16;



	/*glBegin(GL_LINE_STRIP);

	Tool::colors[0].Dim(.2)();

	glVertex3f(this->rect.x,-this->rect.y-this->rect.height,.5);
	glVertex3f(this->rect.x,-this->rect.y,.5);
	glVertex3f(this->rect.x+this->rect.width,-this->rect.y,.5);
	glVertex3f(this->rect.x+this->rect.width,-this->rect.y-this->rect.height,.5);
	glVertex3f(this->rect.x,-this->rect.y-this->rect.height,.5);
	glEnd();
	*/

	int hilite = this->rect.contains(mouse);


	this->rect.border = Tool::INSET | (hilite?Tool::FILL:0);
	this->rect.Draw(0);


	glBegin(GL_QUAD_STRIP);
	// if (col)
	// glColor4f(1., 1., 1., 1.);
	// else
	Tool::colors[0].Dim(.2)();

	int xx = this->rect.x+this->rect.width-icon;
	//int yy = -(this->rect.y);
	//int h = this->rect.height-1-10;

	for (int i=-1;i<=6;i++) {
		color4 co(0,0,1,1);

		int ii = i;
		if (ii<0) ii=0;
		if (ii>5) ii=5;
		if (this->ref==0) {
			co.SetHue(1-ii/5.);
			co.Alpha(.7)();
		}
		else if (this->ref==1) {
			co.SetGrey(1-ii/5.);
			co.Alpha(.7)();
		}
		else
			co.Alpha(1-ii/5.)();


		int ht = getht(i/5.);
		if (i==-1) ht = -this->rect.y-this->rect.padding+2;
		if (i==6) ht = -this->rect.y-(this->rect.height)+this->rect.padding-1;


		//glVertex3f(xx+icon,yy-h*i/5.-5,.5);
		//glVertex3f(xx,yy-h*i/5.-5,.5);

		glVertex3f(xx+icon,ht,.5);
		glVertex3f(xx,ht,.5);
	}
	glEnd();






	this->rect.border =0;
	this->rect.Draw();

	/*double mm = (int)(minval*100+.5)*.01;
	std::string t = ToString(mm);
	int tl = Tool::font->StringLength(t.c_str());
	glColor4f(0.,0.,0.,.7);
	Tool::font->printString(t.c_str(), this->rect.x+ this->rect.width - icon - tl-2,
	      -this->rect.y-this->rect.height+this->rect.padding,.4);

	mm = (int)(maxval*100+.5)*.01;
	t = ToString(mm);
	tl = Tool::font->StringLength(t.c_str());
	glColor4f(0.,0.,0.,.7);
	Tool::font->printString(t.c_str(), this->rect.x+ this->rect.width - icon - tl-2,
	      -this->rect.y-12+this->rect.padding,.4);
	*/

int sstep = step;
	if (sstep<1) sstep = 1;
	
	if (this->rect.height/sstep<15) sstep = (this->rect.height)/15;
	if (sstep<1) sstep = 1;
	if (sstep>100) sstep = 100;

	for (int i=0;i<=sstep;i++) {

		double mm = ((maxval-minval)*i/(double) sstep);

		double vv = unitscale ? ((mm+minval)*unitscale) : (mm+minval) ;
		
		vv = (int)(vv*100+.5)*.01;
		
		std::string t = ToString( vv );
		int  tl = Tool::font->StringLength(t.c_str());
		glColor4f(0.,0.,0.,.7);
		int ht = getht((sstep-i)/(double)sstep);

		Tool::font->printString(t.c_str(), xx-tl-2, ht-5, .6);


		glBegin(GL_LINE_STRIP);
		glVertex3f(xx,ht,.6);
		glVertex3f(xx+3,ht,.6);
		glEnd();

	}

}

template <>
int ToolScale<int>::ToValueType(const std::string &c) const {
	return ToInt(c);
}



template class ToolScale<int>
;
