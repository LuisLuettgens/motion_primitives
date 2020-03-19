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
#include "toolbutton.h"

#include "GL/gl.h"
#include "conversion.h"
#include "../base/point.h"


ToolStateButton::ToolStateButton(int &s, int id_, int image_, int x_, int y_, int modes_)
		: ToolButton(id_,image_,x_,y_), state(s), modes(modes_) {



	//modes = 0; // also nur click ohne state
	//modes = 1; // an + aus
	//modes >= 2; // Bilder durchwahl
}


ToolButton::ToolButton(int id_, int image_, int x_, int y_) {

	imagex = image_%8;
	imagey = image_/8;
	id = id_;
	x= x_;
	y=y_;

	width=20;
	height=20;
	col=0;

	//modes = 0; // also nur click ohne state
	//modes = 1; // an + aus
	//modes >= 2; // Bilder durchwahl
}

int ToolButton::imagesel() {
	return 0;
}

void ToolButton::Draw(const Point<int> &mouse) {
	int icon = 16;
	int s = imagesel();
	int p = pressed();
	

	float tx1 = (imagex+s)/8.;
	float ty1 = imagey/8.;
	float tx2 = (imagex+s+1)/8.;
	float ty2 = (imagey+1)/8.;


	/*glBegin(GL_LINE_STRIP);

	Tool::colors[0].Dim(.2)();

	glVertex3f(x,y+height,.5);
	glVertex3f(x,y,.5);
	glVertex3f(x+width,y,.5);
	glVertex3f(x+width,y+height,.5);
	glEnd();*/

	int hilite = contains(mouse.x,mouse.y);

	Tool::drawRect(Point<int>(x,-y-height+1), Point<int>(x+width-1,-y),.3,
	               ( hilite?Tool::FILL:Tool::FILL2) | (p?Tool::INSET:Tool::OUTSET) );


	/*colors[1]();
		glBegin(GL_LINE_STRIP);
		glVertex3f(p2.x+1, -p1.y, z+.1);
		glVertex3f(x, -y, z+.1);
		glVertex3f(x, -p2.y, z+.1);
		glEnd();

		colors[1]();
		glBegin(GL_LINE_STRIP);
		glVertex3f(x, -p2.y, z+.1);
		glVertex3f(p2.x, -p2.y, z+.1);
		glVertex3f(p2.x, -p1.y, z+.1);
		glEnd();*/

	glEnable(GL_TEXTURE_2D);

	glBegin(GL_QUADS);
	if (col)
		glColor4f(1., 1., 1., 1.);
	else
		Tool::colors[0].Dim(dim())();

	glTexCoord2f(tx1,ty1);
	glVertex3f(x+2,y+icon+2,.5);
	glTexCoord2f(tx1,ty2);
	glVertex3f(x+2,y+2,.5);
	glTexCoord2f(tx2,ty2);
	glVertex3f(x+2+icon,y+2,.5);
	glTexCoord2f(tx2,ty1);
	glVertex3f(x+2+icon,y+2+icon,.5);
	glEnd();
	glDisable(GL_TEXTURE_2D);


}


bool ToolButton::contains(int x_,int y_) {

	if (x_>x && x_<x+width && y_>y && y_<y+height)
		return true;
	return false;
}

float ToolStateButton::dim() {
	if (modes>1)
		return .4;
	return state?.5:.1;
}

float ToolButton::dim() {
	return .4;
}

int ToolStateButton::imagesel() {
	if (modes==1)
		return 0;
	return state;
}

int ToolStateButton::pressed() {
	if (modes==1)
		
	return state;
	return 0;
}
