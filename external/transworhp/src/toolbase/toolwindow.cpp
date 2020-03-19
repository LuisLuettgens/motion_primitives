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

#include "toolwindow.h"

#include "GL/gl.h"
#include "conversion.h"

using namespace std;

namespace tw {

int ToolWindow::scrollht = 16;
int ToolWindow::scrollwd = 7;


ToolWindow::ToolWindow(const std::string &t) : title(t) {

	moving = 0;
	width = 150+10;
	height = 140;

	x = 0;
	y = 0;

	open = 1;
	theight = 18;

	closeme = false;
	align=0;

	visible = 1;


	setActiveModeTo = -1;

	scrollmax = height;
	scrollcur = 0;
	scrolling = false;

	resizeable = 0;
	resizing = false;
	minsizex = width;
	minsizey = height;
}


ToolWindow::~ToolWindow() {
	// clearpointer(buttons);

	clearpointer(icons);

}

int ToolWindow::posleft = 21;
int ToolWindow::posright = 21;

void ToolWindow::Init() {}

void ToolWindow::InitPos(int align_, int open_) {

	open = open_;
	align = align_;

	//std::cout << align << open << " " << ww << " " <<x << y << title << endl;
}

void ToolWindow::move(int ww, int /*hh*/) {

	if (align==1) {
		x = 0;
		y = posleft;
		posleft += (open? height:theight)+2;

	}
	if (align==2) {
		x = ww - width -1;
		y = posright;
		posright += (open? height:theight)+2;
	}

}



ToolIcon::ToolIcon(ToolIcon::action id_, int image_) {

	tx = image_%8;
	ty = image_/8;

	id = id_;

}


void ToolIcon::Draw(int width, int pos, bool moving, bool flag) {

	int a = 0;
	if (flag && id==OPEN)
		a++;
	int icon = 16;
	float tx1 = (tx+a)/8.f;
	float ty1 = ty/8.f;
	float tx2 = (tx+1+a)/8.f;
	float ty2 = (ty+1)/8.f;
	glEnable(GL_TEXTURE_2D);
	Tool::texture->Bind();
	glBegin(GL_QUADS);

	glColor4fv(Tool::colors[0].Dim(moving?.6f:.4f).GetData());

	glTexCoord2f(tx1,ty1);
	glVertex3f((float)(width+pos*icon),icon,.5f);
	glTexCoord2f(tx1,ty2);
	glVertex3f(width+pos*icon,0,.5f);
	glTexCoord2f(tx2,ty2);
	glVertex3f(width+(pos+1)*icon,0,.5f);
	glTexCoord2f(tx2,ty1);
	glVertex3f(width+(pos+1)*icon,icon,.5f);
	glEnd();
	glDisable(GL_TEXTURE_2D);
}


void ToolWindow::AddIcon(ToolIcon::action id) {

	int image = 0;
	if (id==ToolIcon::OPEN)
		image = 0;
	if (id==ToolIcon::CLOSE)
		image = 7;
	if (id==ToolIcon::WRITEPNG)
		image = 3+8*3;
	if (id==ToolIcon::WRITETXT)
		image = 4+8*3;

	icons.push_back(new ToolIcon(id,image));

}


void ToolWindow::IconAction(ToolIcon::action id) {

	if (id == ToolIcon::OPEN) {
		open ^=1;
	}

	if (id == ToolIcon::WRITETXT) {
		WriteText();
	}

	if (id == ToolIcon::WRITEPNG) {
		WriteImage();
	}

	if (id == ToolIcon::CLOSE) {
		closeme = true;
	}

}




void ToolWindow::Draw(bool isActive, int /*HT*/, const Point<int> &/*mouse*/, bool border) {

	if (!visible)
		return;

	glDisable(GL_LINE_SMOOTH);

	// Tiefe von 0 bis 1

	// Titelzeile
	if (border) {

		glBegin(GL_QUADS);

		//Tool::colors[0].Dim(moving?.1:.0)();

		//glVertex3f(x+2,-y,0);
		//glVertex3f(x+2,-y-theight,0);

		//Tool::colors[0].Dim(moving?.15:.05)();

		//glVertex3f(x+width-2,-y-theight,0);
		//glVertex3f(x+width-2,-y,0);


		// Textfeld

		if (open) {
			glColor4fv(Tool::colors[5].Dim(.0f).GetData());

			glVertex3f(x,-y,.0f);
			glVertex3f(x,-y-height,.0f);

			glColor4fv(Tool::colors[5].Dim(.05f).GetData());

			glVertex3f(x+width,-y-height,.0f);
			glVertex3f(x+width,-y,.0f);

		}
		else  {
			glColor4fv(Tool::colors[5].Dim(.0f).GetData());

			glVertex3f(x,-y,.0f);
			glVertex3f(x,-y-theight,.0f);

			glColor4fv(Tool::colors[5].Dim(.05f).GetData());

			glVertex3f(x+width,-y-theight,.0f);
			glVertex3f(x+width,-y,.0f);

		}
		glEnd();

	}

	Tool::drawRect(Point<int>(x+2,y+2), Point<int>(x+width-2,y+theight-2), .1, (moving?Tool::FILLDARK : Tool::FILL2) | Tool::INSET);
	// Tool::drawRect(Point(x+2,y+2), Point(x+width-2,y+theight-2), .1, (isActive?Tool::FILLHI : Tool::FILL2) | Tool::INSET);


	/*if (isActive)
	 glLineWidth(2);
	else {
	 glLineWidth(1);
	}*/
	glLineWidth(1);

	glDisable(GL_LINE_SMOOTH);

	if (border) {

		int h = height;
		if (!open)
			h=theight;

		Tool::drawRect(Point<int>(x,y),Point<int>(x+width,y+h),0., Tool::OUTSET);

		if (isActive) {
			Tool::drawRect(Point<int>(x-1,y-1),Point<int>(x+width+1,y+h+1),0., Tool::SOLID);
		}

		/*
		Tool::colors[1]();
		 glBegin(GL_LINE_STRIP);
		 glVertex3f(x,-y,.1);
		 glVertex3f(x+width+.1,-y,.1);
		Tool::colors[3]();
		 if (open) {
		  glVertex3f(x+width+.1,-y-height,.1);
		  glVertex3f(x,-y-height,.1);
		 } else {
		  glVertex3f(x+width+.1,-y-theight,.1);
		  glVertex3f(x,-y-theight,.1);
		 }
		 glVertex3f(x,-y,.1);
		 glEnd();
		}*/
	}


	glTranslatef(x,-y-theight+1,0);


	int w = Tool::font->StringLength(title.c_str());

	if (isActive) {

		glColor4fv(Tool::colors[3].GetData());
		Tool::font->printString(title.c_str(), (width-w)/2-1, 3+1, .4);
	}
	glColor4fv(Tool::colors[6].Alpha(.8).GetData());
	Tool::font->printString(title.c_str(), (width-w)/2, 3, .5);

	glLineWidth(1);

	/* if (scrollmax > height && open) {
	  double ht = (double)(height-scrollht-theight)/(scrollmax - height)*scrollcur;//+theight;
	 
	  Tool::colors[0].alpha(.5);
	 
	  glBegin(GL_QUADS);
	 
	  glVertex3f(width-scrollwd,-height-theight,.05);
	  glVertex3f(width+.1,-height-theight,.05);
	  Tool::colors[0].alpha(.2);
	  glVertex3f(width+.1,0,.05);
	  glVertex3f(width-scrollwd,0,.05);
	  glEnd();
	 
	 
	  int icon = 16;
	  int tx = 7,ty = 2;
	  float tx1 = tx/8.;
	  float ty1 = ty/8.;
	  float tx2 = (tx+1)/8.;
	  float ty2 = (ty+1)/8.;
	  glEnable(GL_TEXTURE_2D);
	  Tool::texture->Bind();
	  glBegin(GL_QUADS);
	  //glBegin(GL_LINE_STRIP);
	 
	  Tool::colors[0].dim(moving?.6:.4);
	 
	  glTexCoord2f(tx1,ty1);
	  glVertex3f(width-icon,-ht,.1);
	  glTexCoord2f(tx1,ty2);
	  glVertex3f(width-icon,-ht-icon,.1);
	  glTexCoord2f(tx2,ty2);
	  glVertex3f(width,-ht-icon,.1);
	  glTexCoord2f(tx2,ty1);
	  glVertex3f(width,-ht,.1);
	  glEnd();
	  glDisable(GL_TEXTURE_2D);
	 }*/


	if (resizeable && open && isActive) {

		glColor4fv(Tool::colors[1].GetData());
		glBegin(GL_LINES);

		int x = width;
		int y = height-theight;

		for (int i=2;i<11;i+=3) {
			glVertex3f(x-i, -y, .2);
			glVertex3f(x, -y+i, .2);
		}
		glEnd();

		/*	Tool::drawRect(Point<int>(width-scrollwd,height-theight-scrollwd),Point<int>(width-2,height-theight-1),0.1, Tool::INSET);
		*/

	}


	if (scrollmax > height && open) {
		int ht = (int)((double)(height-scrollht-theight-3)/(scrollmax - height)*scrollcur)+1;//+theight;

		Tool::drawRect(Point<int>(width-scrollwd,0),Point<int>(width-2,height-theight-1),0.1, Tool::INSET);


		Tool::drawRect(Point<int>(width-scrollwd+1,ht),Point<int>(width-3,ht+scrollht),0.1, Tool::OUTSET | Tool::FILLDARK);



		/*
		  glBegin(GL_QUAD_STRIP);

		  Tool::colors[0].Alpha(.2)();
		  glVertex3f(width,0,.05);
		  glVertex3f(width-scrollwd,0,.05);

		  Tool::colors[0].Alpha(.2)();
		  glVertex3f(width,-ht,.1);
		  glVertex3f(width-scrollwd,-ht,.1);


		  Tool::colors[0].Alpha(1.)();
		  glVertex3f(width,-ht-scrollht/2,.1);

		  glVertex3f(width-scrollwd,-ht-scrollht/2,.1);

		  Tool::colors[0].Alpha(.5)();
		  glVertex3f(width,-ht-scrollht,.1);
		  glVertex3f(width-scrollwd,-ht-scrollht,.1);

		  Tool::colors[0].Alpha(.5)();
		  glVertex3f(width,-height+theight,.05);
		  glVertex3f(width-scrollwd,-height+theight,.05);


		  glEnd();
		*/

		/* int icon = 16;
		 int tx = 7,ty = 2;
		 float tx1 = tx/8.;
		 float ty1 = ty/8.;
		 float tx2 = (tx+1)/8.;
		 float ty2 = (ty+1)/8.;
		 glEnable(GL_TEXTURE_2D);
		 Tool::texture->Bind();
		 glBegin(GL_QUADS);
		 //glBegin(GL_LINE_STRIP);

		 Tool::colors[0].dim(moving?.6:.4);

		 glTexCoord2f(tx1,ty1);
		 glVertex3f(width-icon,-ht,.1);
		 glTexCoord2f(tx1,ty2);
		 glVertex3f(width-icon,-ht-icon,.1);
		 glTexCoord2f(tx2,ty2);
		 glVertex3f(width,-ht-icon,.1);
		 glTexCoord2f(tx2,ty1);
		 glVertex3f(width,-ht,.1);
		 glEnd();
		 glDisable(GL_TEXTURE_2D);*/
	}


	if (border) {
		std::vector<ToolIcon*>::iterator it = icons.begin();
		int ii=0;
		for (;it!=icons.end();it++) {
			(*it)->Draw(width,--ii,moving,open);
		}
	}
}

std::string ToolWindow::GetFilename() {
	return "test";
}


void ToolWindow::WriteText() {}
void ToolWindow::WriteImage() {}


bool ToolWindow::KeyboardFunc(SDL_Keysym &/*Keysym*/) {
	return true;
}


void ToolWindow::Resize() {}

bool ToolWindow::mousebutton(int /*mx*/, int /*my*/, int /*button*/, int /*type*/) {
	return true;
}

bool ToolWindow::mousemotion(int /*mx*/, int /*my*/) {
	return true;
}


bool ToolWindow::MouseMotion(bool isActive, SDL_MouseMotionEvent &m) {

	if (!visible)
		return false;

	if (isActive && scrolling) {
		if (open) {
			//if (m.x>x && m.x<x+width && m.y>y+theight && m.y<y+height && open) {
			if (m.state & SDL_BUTTON(1)) { // && m.x>width-scrollwd) {

				scrollcur = (int)((m.y-y-scrollht/2 -theight) * (scrollmax-height)
								  / (double)(height-scrollht-theight));
				if (scrollcur<0)
					scrollcur = 0;
				if (scrollcur>scrollmax-height)
					scrollcur = scrollmax-height;
			}
		}
	}


	if (isActive && resizing) {
		if (open) {
			//if (m.x>x && m.x<x+width && m.y>y+theight && m.y<y+height && open) {
			if (m.state & SDL_BUTTON(1)) { // && m.x>width-scrollwd) {

				if (resizeable & 1) {
					height = (int)(m.y-y);
if (height<minsizey) height =minsizey;
					scrollcur = 0;
					scrollmax = height;
				}
				if (resizeable & 2) {
					width = (int)(m.x-x);
if (width<minsizex) width = minsizex;
				}

				Resize();
				/*if (scrollcur<0)
				 scrollcur = 0;
				if (scrollcur>scrollmax-height)
				 scrollcur = scrollmax-height;
				*/
			}
		}
	}


	if (moving) {
		int tmpx = m.x-mouse.x;
		int tmpy = m.y-mouse.y;

		x = tmpx;
		y = tmpy;
		return true;
	}

	if (m.x>x && m.x<x+width && m.y>y+theight && m.y<y+height && open) {

		int mx = m.x-x;
		int my = m.y-y-theight;

		mousemotion(mx,my+scrollcur);
	}

	return false;
}


bool ToolWindow::MouseButton(SDL_MouseButtonEvent &m) {

	if (!visible)
		return false;


	if ((m.button == 1 || m.button==4 || m.button==5) && m.type == SDL_MOUSEBUTTONUP) {
		if (moving)
			moving=0;

		scrolling = false;
		resizing = false;
	}


	// Scrollen per Mausrad
	if ((m.button == 4 || m.button == 5) && m.type == SDL_MOUSEBUTTONDOWN && scrollmax>height) {

		if (m.x>x+width-10 && m.x<x+width && m.y>y+theight && m.y<y+height && open) {

			if (m.button == 5)
				scrollcur+=15;
			if (m.button == 4)
				scrollcur-=15;
			if (scrollcur>scrollmax-height)
				scrollcur = scrollmax-height;
			if (scrollcur<0)
				scrollcur = 0;

			return true;
		}
	}

	if (resizeable) {
// Scrollen per Leistenklick
		if (m.button == 1 && m.type == SDL_MOUSEBUTTONDOWN) {
			if (m.x>x+width-10 && m.x<x+width && m.y>y+height-10 && m.y<y+height && open) {
				resizing = true;
			}
		}
	}


	// Scrollen per Leistenklick
	if (m.button == 1 && m.type == SDL_MOUSEBUTTONDOWN  && scrollmax>height) {
		if (m.x>x+width-scrollwd && m.x<x+width && m.y>y+theight && m.y<y+height && open) {

			scrollcur = (int)((m.y-y-scrollht/2 -theight) * (scrollmax-height)
							  / (double)(height-scrollht-theight));
			if (scrollcur>scrollmax-height)
				scrollcur = scrollmax-height;
			if (scrollcur<0)
				scrollcur = 0;

			scrolling = true;
		}
	}


	// Rechte Maustaste weitergeben
	if ((m.button == 3) && m.type == SDL_MOUSEBUTTONDOWN) {

		if (m.x>x && m.x<x+width && m.y>y+theight && m.y<y+height && open) {

			int mx = m.x-x;
			int my = m.y-y-theight;

			mousebutton(mx,my+scrollcur,m.button,m.type);

			return true;
		}
	}

	// Linke Maustaste
	if ((m.button == 1 || m.button==4 || m.button==5) && m.type == SDL_MOUSEBUTTONDOWN) {

		// in Fenster
		if (m.x>x && m.x<x+width && m.y>y && m.y<y+theight) {

			if (setActiveModeTo>=0)
				Tool::activeMode = setActiveModeTo;

			int ic =0;
			std::vector<ToolIcon *>::iterator it = icons.begin();
			for (;it!=icons.end();it++,ic++) {

				if (m.x>x+width-16-16*ic && m.x<x+width-16*ic
						&& m.y>y && m.y<y+16) {

					IconAction((*it)->id);
					return true;
				}
			}

			if (m.y>y && m.y<y+theight) {
				mouse = Point<int>(m.x-x, m.y-y);

				moving = 1;
				return true;
			}

			return true;

		}
		else if (m.x>x && m.x<x+width && m.y>y+theight && m.y<y+height && open) {

			if (setActiveModeTo>=0)
				Tool::activeMode = setActiveModeTo;

			int mx = m.x-x;
			int my = m.y-y-theight;

			mousebutton(mx,my+scrollcur,m.button,m.type);
			/*

			   std::vector<ToolButton*>::const_iterator it = buttons.begin();
			   for (;it!=buttons.end();it++) {

			    if ((*it)->contains(mx,my-theight)) {
			     gl->Message((*it)->id);
			     return true;
			    }
			   }

			  // cout << mx << " " << my << endl;

			   std::vector<ToolItem*>::iterator it2 = items.begin();
			   for (;it2!=items.end();it2++) {
			    if ((*it2)->TakesKey()) {
			     if ((*it2)->contains(mx,my)) {
			      (*it2)->MouseClick(m.button,m.type,Point(mx,my));
			      activeItem = (it2);
			      return true;
			     }
			    }
			   }

			*/
			return true;
		}
		else {
			moving = 0;
		}
	}



	return false;
}



/*
bool ToolWindow::SpecialKeys(Uint8 *keys) {

 if (activeItem) {
  activeItem->SpecialKey(keys);
 }

 return true;
}*/

void ToolWindow::innerProj(int HT) {

	glViewport(x, HT-y-height+1, width, height-theight-1);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(x, x+width, -y-height+1,-y-theight, -12.0, 2.0);
	glMatrixMode(GL_MODELVIEW);

	glTranslatef(0,scrollcur,0);
}












void ToolWindow::RenderWindow() {

	glClearColor(1.f, 1.f, 1.f, 1.f);

	// Clear the window with current clearing color
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	// glLightfv(GL_LIGHT0, GL_POSITION, light_pos);

	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);

	glEnable(GL_ALPHA_TEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glClear(GL_DEPTH_BUFFER_BIT);

	//displayText3d(screen,xtmp, ytmp);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glDisable(GL_LIGHTING);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0, width, -height, 0, -12.0, 2.0);
	glMatrixMode(GL_MODELVIEW);

	int x0 = x, y0 = y;
	x=0;
	y=0;
	glPushMatrix();
	//glTranslatef(-x,y,0);
	Draw(false,height,mouse,false);
	glPopMatrix();

	x=x0;
	y=y0;

	glEnable(GL_LIGHTING);

	//displayText(screen);

	glDisable(GL_ALPHA_TEST);
	glDisable(GL_BLEND);

}

}
