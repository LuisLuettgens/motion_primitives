#include "viewport.h"

namespace tw {

void Viewport::SetLeft() {
	glViewport(r.x/2,r.y,r.width/2,r.height);
}

void Viewport::SetRight(int w) {
	glViewport(r.x/2+w,r.y,r.width/2,r.height);
}

void Viewport::Set() {
	glViewport(r.x,r.y,r.width,r.height);
}

void Viewport::Position() {
	glTranslatef(0.f,0.f,-(float)tau);
	glRotatef(-(float)yang,1.f,0.f,0.f);
	glRotatef(-(float)xang,0.f,1.f,0.f);
glTranslatef(-(float)dx,-(float)dy,0.f);
	glRotatef(-(float)roll,0.f,0.f,1.f);
glRotatef(-(float)neigung,1.f,0.f,0.f);

}

/*
void Viewport::Frame(int h) {
	glVertex3f(x,       -h+y,       0);
	glVertex3f(x,       -h+y+height,0);
	glVertex3f(x+width, -h+y+height,0);
	glVertex3f(x+width, -h+y,       0);
	glVertex3f(x,       -h+y,       0);
}*/

}
