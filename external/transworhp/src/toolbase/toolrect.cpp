#include "toolrect.h"

#include <iostream>

using namespace std;

namespace tw {

bool ToolRect::contains(int x_,int y_) const {

	if (x_>x && x_<x+width && y_>y && y_<y+height)
		return true;

	return false;
}

bool ToolRect::contains(const Point<int> &p) const {

	if (p.x>x && p.x<x+width && p.y>y && p.y<y+height)
		return true;

	return false;
}


void ToolRect::Draw(double b, float addz) const {

	// hmm...
	if (b) {
		glColor4fv(Tool::colors[5].Dim(b).GetData());
		glBegin(GL_QUADS);
		glVertex3f(x,-y-height,.25+addz);
		glVertex3f(x+width,-y-height,.25+addz);
		glVertex3f(x+width,-y,.25+addz);
		glVertex3f(x,-y,.25+addz);
		glEnd();
	}

	if (border) {
		//Tool::colors[5].Dim(.6)();

		Tool::drawRect(Point<int>(x,y),Point<int>(x+width,y+height),.3+addz,border);
		/*
		glBegin(GL_LINE_STRIP);
		glVertex3f(x,-y,.3+addz);
		glVertex3f(x+width,-y,.3+addz);
		glVertex3f(x+width,-y-height,.3+addz);
		glVertex3f(x,-y-height,.3+addz);
		glVertex3f(x,-y,.3+addz);
		glEnd();*/
	}

	glColor4fv(Tool::colors[5].Dim(.8).GetData());

	int a = text.size();

	if (a>0 && a<1000) {
		if (border&256) {
			int ww = Tool::font->StringLength(text.c_str());
			Tool::font->printString(text.c_str(), x+(width-ww)/2 +1 , -y-height+padding, .4+addz,width);

		} else {
			Tool::font->printString(text.c_str(), x+padding+1, -y-height+padding, .4+addz,width);
		}
	}
	else if (a>=1000){

	cout << "ToolRect::Draw() textsize=" << a << endl;
	//cout << text << endl;
	}

}

}
