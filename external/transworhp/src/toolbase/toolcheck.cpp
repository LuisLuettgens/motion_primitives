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
//#include "toolitem.h"
#include "GL/gl.h"

#include "conversion.h"
#include "toolcheck.h"



template <class T>
std::string ToolCheck<T>::GetString() {
	return ToString(this->value);
}

template<>
void ToolCheck<int>::SetString(const std::string &s) {
	value = ToInt(s);
}


template <class T>
void ToolCheck<T>::MouseClick(int button, int state, const Point<int> &p) {

	if (button==1 && state==SDL_MOUSEBUTTONDOWN) {

		this->value = 1 - this->value;

		this->Apply();
	}
}




template <class T>
void ToolCheck<T>::Key(SDL_Keysym &keysym) {

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
void ToolCheck<T>::Draw(bool isActive, const Point<int> &mouse) {

	int hilite = this->rect.contains(mouse);

	int icon = 14;
	int s = this->value?1:0;

	float tx1 = (s)/8.+1/(8.*icon);
	float ty1 = 3/8.+1/(8.*icon);
	float tx2 = (s+1)/8.-1/(8.*icon);
	float ty2 = (3+1)/8.-1/(8.*icon);


	ToolRect r2(this->rect);
	int xx = this->rect.x+this->rect.width-icon;
	int yy = -(this->rect.y);
	
	r2.border = Tool::INSET | (hilite?Tool::FILL:0);
	r2.x =  this->rect.x+this->rect.width-icon;
	r2.width = icon;
	r2.text = "";
	r2.Draw();
	//Tool::drawRect(Point(xx,-yy),Point(xx+icon,-yy+icon),.4,Tool::INSET | (hilite?Tool::FILL:0));


	//Tool::colors[1]();
	
	glEnable(GL_TEXTURE_2D);

	Tool::texture->Bind();


	glBegin(GL_QUADS);
	//	if (col)
	//	glColor4f(1., 1., 1., 1.);
	//	else


	Tool::colors[5].Dim(.8)();

//	Tool::colors[0].Dim(.7)();
	glTexCoord2f(tx2,ty2);
	glVertex3f(xx+icon,yy-icon,.5);
	glTexCoord2f(tx2,ty1);
	glVertex3f(xx+icon,yy,.5);
	glTexCoord2f(tx1,ty1);
	glVertex3f(xx,yy,.5);
	glTexCoord2f(tx1,ty2);
	glVertex3f(xx,yy-icon,.5);
	glEnd();
	glDisable(GL_TEXTURE_2D);


	this->rect.border = 0;
	this->rect.width-=20;
	this->rect.Draw(0);
	this->rect.width+=20;



}

template <>
int ToolCheck<int>::ToValueType(const std::string &c) const {
	return ToInt(c);
}



template class ToolCheck<int>
;
