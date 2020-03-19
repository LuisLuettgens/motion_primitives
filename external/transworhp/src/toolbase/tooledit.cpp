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

#include "../base/defines.h"
#include "conversion.h"
#include "tooledit.h"





template <class T>
std::string ToolEdit<T>::GetString() {
	return ToString(this->value);
}


template<>
std::string ToolEdit<double>::GetString() {

	if (scale)
		return ToString(this->value* *scale);
	else
		return ToString(this->value);
}

template<>
std::string ToolEdit<std::string>::GetString() {
	return value;
}

template<>
void ToolEdit<int>::SetString(const std::string &s) {
	value = ToInt(s);
}

template<>
void ToolEdit<double>::SetString(const std::string &s) {
	value = ToDouble(s);
}

template<>
void ToolEdit<std::string>::SetString(const std::string &s) {
	value = s;
}


template <class T>
void ToolEdit<T>::Discard() {
	ToolData<T>::Discard();
	edittext = GetString();
}


/** Spezialisierungen von AnyEdit */
template <class T>
void ToolEdit<T>::MouseClick(int button, int state, const Point<int> &p) {

	if (button==1 && state==SDL_MOUSEBUTTONDOWN) {
		//  edittext = GetString();
		//  cursor=edittext.size();
		//  editmode=1;
		// }

		for (unsigned int i=0;i<edittext.size();i++) {
			std::string a(edittext,0,i);

			int aa=Tool::font->StringLength(a.c_str());
			int bb=p.x-this->rect.x-3-this->rect.padding;

			if (aa>bb) {
				cursor=i;
				return;
			}
		}

		cursor=edittext.size();

	}
}


template <>
bool ToolEdit<std::string>::CheckKey(SDL_Keycode &c) {

	if (c>' ')
		return true;
	return false;
}

template <>
bool ToolEdit<int>::CheckKey(SDL_Keycode &c) {

	if (c>='0' && c<='9')
		return true;
	if (c>=SDLK_KP_0 && c<=SDLK_KP_9) {
		c = (SDL_Keycode)(c- (SDLK_KP_0 - '0'));
		return true;
	}
	if (c=='-' || c=='+')
		return true;
	return false;
}

template <>
bool ToolEdit<double>::CheckKey(SDL_Keycode &c) {

	if (c>='0' && c<='9')
		return true;
	if (c>=SDLK_KP_0 && c<=SDLK_KP_9) {
		c = (SDL_Keycode)(c- (SDLK_KP_0 - '0'));
		return true;
	}
	if (c=='-' || c=='+' || c=='.')
		return true;
	return false;
}


template <class T>
void ToolEdit<T>::Key(SDL_Keysym &keysym) {

	SDL_Keycode key = keysym.sym;
	// SDLMod mod = keysym.mod;

	//std::cout << (int)c << std::endl;

	switch (key) {
	case 8: // <--
		if (cursor>=1) {
			edittext.erase(cursor-1,1);
			cursor--;
		}
		break;

	case 13: // RETURN
		// editmode = 0;
		this->value = ToValueType(edittext);
		break;

	case 27: // ESC
		// editmode = 0;
		break;

	case 127: // ENTF
		if (cursor>=0 && cursor<(int)edittext.size()) {
			edittext.erase(cursor,1);
			//std::cout << (int)key << std::endl;
		}
		break;

	case SDLK_LEFT:
		cursor--;
		if (cursor<0)
			cursor=0;
		break;

	case SDLK_RIGHT:
		cursor++;
		if (cursor>(int)edittext.size())
			cursor=edittext.size();
		break;
	default:
		if (CheckKey(key)) {
			edittext.insert(cursor,1,key);
			//std::cout << "key" << key << std::endl;
			this->value = ToValueType(edittext);
			cursor++;
		}
	}

	if (applyimmediately) {
		this->value = ToValueType(edittext);

		this->Apply();
	}
}




template <class T>
void ToolEdit<T>::Draw(bool isActive, const Point<int> &mouse) {



	/* this->back.z = 0.1f;
	 this->back.SetCol(2);
	 this->back.border = 0;
	 this->back.mode = 0;
	 this->back.centertext=0;
	 //std::cout << "..." << p << "..." << back.t << " "<< back.b<<" " << back.l << " "
	 //<<  back.r << std::endl;
	 *if (this->back.Contains(p)) {
	  this->back.mode =1;
	 }
	 else {
	  if (editmode==1) {
	   editmode = 0;
	   this->value = ToValueType(edittext);
	  }
	 }*/

	// if (editmode) {
	//this->back.SetCol(1);

	int hilite = this->rect.contains(mouse);

	this->rect.border = Tool::INSET | (hilite?Tool::FILL:0);

	this->rect.text = edittext;
	// }
	// else {
	//  this->rect.text = GetString();
	// }

	this->rect.Draw();

	if (isActive) {
		int i = Tool::font->StringLength(std::string(edittext,0,cursor).c_str())
		        ;

		/* glColor4f(1.,0.,0.,1.);
		 glRasterPos3f(20,-20,back.z+2.);
		 glDialog::rf->printString((ToString(cursor) + " " +
		  ToString(i)).c_str());
		*/

		glColor4f(1.,0.,0.,1.);
		//glRasterPos3f(back.l+2+i,back.b+4-1,back.z+4.);
		Tool::font->printString("|",
		                        this->rect.x+i+this->rect.padding,
		                        -this->rect.y-this->rect.height+2+this->rect.padding,.6);
	}

}

template <>
int ToolEdit<int>::ToValueType(const std::string &c) const {
	return ToInt(c);
}

template <>
double ToolEdit<double>::ToValueType(const std::string &c) const {
	if (scale)
		if (*scale)
			return ToDouble(c)*(1. / *scale);

	return ToDouble(c);
}

template <>
std::string ToolEdit<std::string>::ToValueType(const std::string &c) const {
	return c;
}



template class ToolEdit<double>
;
template class ToolEdit<std::string>
;
template class ToolEdit<int>
;


