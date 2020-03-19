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
#include "toolslider.h"



template <class T>
std::string ToolSlider<T>::GetString() {
	return ToString(this->value);
}

template <>
std::string ToolSlider<double>::GetString() {
	return ToString(this->value);
}

template<>
void ToolSlider<int>::SetString(const std::string &s) {
	value = ToInt(s);
}

template<>
void ToolSlider<double>::SetString(const std::string &s) {
	value = ToDouble(s);
}

template <class T>
std::string ToolSlider<T>::GetRefString() {
	return ToString(this->ref);
}

template <>
std::string ToolSlider<double>::GetRefString() {
	//return ToString(this->ref,0,precision);
	/*char buf[40];
	char buf2[40];
	sprintf(buf,"%%.%df",precision);
	double a = this->ref;
	sprintf(buf2,buf,a);
	return std::string(buf2);*/
	return ToString(this->ref,0,precision);
}






/** Spezialisierungen von AnyEdit */
template <class T>
void ToolSlider<T>::MouseClick(int button, int state, const Point<int> &p) {
/*
	if (button==1 && state==SDL_MOUSEBUTTONDOWN) {
		//  edittext = GetString();
		//  cursor=edittext.size();
		//  editmode=1;
		// }

	}*/
}



template <>
bool ToolSlider<int>::CheckKey(SDL_Keycode &c) {

	if (c>='0' && c<='9')
		return true;
	if (c=='-' || c=='+')
		return true;
	return false;
}

template <>
bool ToolSlider<double>::CheckKey(SDL_Keycode &c) {

	if (c>='0' && c<='9')
		return true;
	if (c=='-' || c=='+' || c=='.')
		return true;
	return false;
}


template <class T>
void ToolSlider<T>::Key(SDL_Keysym &keysym) {
/*
	SDLKey key = keysym.sym;
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
			std::cout << (int)key << std::endl;
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
			this->value = ToValueType(edittext);
			cursor++;
		}
	}*/
}




template <class T>
void ToolSlider<T>::Draw(bool isActive, const Point<int> &mouse) {

	this->rect.border = Tool::INSET;
	this->rect.Draw();

	Tool::colors[5].Dim(.20f)();

	int w = 50;

	int a = (int)(this->rect.x + (this->rect.width-w-3)*(this->ref-low)/(high-low));
	//a = 20;
	
	ToolRect z(a+1,this->rect.y + 1,w+1,this->rect.height - 2);
	z.text=GetRefString();
	z.padding=1;
z.border = Tool::OUTSET | Tool::CENTER;
	z.Draw(.3f,.1f);

}

template <>
int ToolSlider<int>::ToValueType(const std::string &c) const {
	return ToInt(c);
}

template <>
double ToolSlider<double>::ToValueType(const std::string &c) const {
	return ToDouble(c);
}


template class ToolSlider<double>;
template class ToolSlider<int>;
