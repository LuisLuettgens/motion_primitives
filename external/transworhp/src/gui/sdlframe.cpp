#include "sdlframe.h"

#include "sdlthread.h"

#include "../base/point.h"

#include "../imaging/imagewriter.h"

#define _USE_MATH_DEFINES
#include <cmath>

#ifndef M_PI
#define M_PI 3.14152926535897932846
#endif

#ifdef WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

using namespace std;

namespace tw {

SDLFrame::SDLFrame(SDLScreen *screen_)
	:	screen(screen_), //fullscreen(false), width(800), height(600),
		//res_x(800), res_y(600), pos_x(0), pos_y(0),
		waittime(0), returnvalue(0),
		imgwriter(nullptr), imgwriterflag(0),
		render(false), busy(false) {

}

SDLFrame::~SDLFrame() {

}


Uint32 sdltimer(Uint32 interval, void *param) {

	SDLFrame *f = (SDLFrame*)param;

	if (f) {
		f->Timer();
	}

	return interval;
}



void SDLFrame::Timer() {

	if (!busy) {
		busy = true;

		const Uint8 *keys = SDL_GetKeyboardState(NULL);
		SpecialKeys(keys);

		TimerFunc();

		//render = true;
		//SDL_SemWait(sem1);


		/*SDL_Event e;
		e.type = SDL_VIDEOEXPOSE;
		SDL_PushEvent(&e);*/

		busy = false;
	}
}


void SDLFrame::MyQuit(int status) {

	returnvalue = status;

	SDL_Event event;
	event.type = SDL_QUIT;
	SDL_PushEvent(&event);
}


void SDLFrame::on_SDL_QUIT() {
	//hier leer
}


int SDLFrame::Loop(int wait, int terminate) {

	waittime = wait;

	//cout << SDL_GetTicks() << "*** SDLFrame::Loop..." << endl;

	/* MR: warum kein cursor anzeigen?
	if (screen->twwindow.fullscreen) {
		SDL_ShowCursor(SDL_DISABLE);
	}*/

	SDL_TimerID my_timer_id = SDL_AddTimer(20, sdltimer, this);
	SDL_Event event;

	bool done = false;
	while (!done) {

		//Idle();

		while (SDL_PollEvent(&event)) {
			switch (event.type) {
			case SDL_WINDOWEVENT:

				switch (event.window.event) {

					case SDL_WINDOWEVENT_RESIZED:
         //   SDL_Log("Window %d resized to %dx%d",
          //          event->window.windowID, event->window.data1,
         //           event->window.data2);

						Reshape(event.window.data1,event.window.data2);
            break;
				//	screen->screen = SDL_SetVideoMode(
				//				 event.resize.w, event.resize.h, 16,
				//				 SDL_OPENGL | SDL_RESIZABLE);
				//	if (screen->screen) {
				//		Reshape(event.resize.w, event.resize.h);
				//	}
				//	else {
				//		// Uh oh, we couldn't set the new video mode?? ;
				//	}
				}

				break;

			case SDL_KEYDOWN:
				KeyboardFunc(event.key.keysym);
				break;

			case SDL_MOUSEBUTTONUP:
			case SDL_MOUSEBUTTONDOWN:
				MouseButton(event.button);
				break;
			
			case SDL_MOUSEWHEEL:
				MouseWheel(event.wheel);
				break;

			case SDL_MOUSEMOTION:
				MouseMotion(event.motion);
				break;

			case SDL_JOYBUTTONUP:
			case SDL_JOYBUTTONDOWN:
				JoystickButtonEvent(event.jbutton);
				break;

			case SDL_JOYAXISMOTION:
				JoystickAxisEvent(event.jaxis);
				break;

			case SDL_JOYHATMOTION:
				JoystickHatEvent(event.jhat);
				break;

			case SDL_QUIT:
				on_SDL_QUIT();
				done = true;
				break;

			case SDL_USEREVENT:{
				Message(reinterpret_cast<intptr_t>(event.user.data1));
			}
				break;

			//case SDL_VIDEOEXPOSE:
			//	break;

			}
		}



		if (screen->thethread && terminate) {
			if (screen->thethread->threadlive>=2) {
				screen->thethread->threadlive++;

				if (screen->thethread->threadlive>=100) {
					//cout << "QUIT THREADLIVE = 1" << endl;
					screen->thethread->threadlive=1;
					done=true;
				}

			}
		}

		if (render) {

		    //	cout << SDL_GetTicks() << "*** SDLFrame::Please Render" << endl;

			render=false;



	RenderScene();


    	//cout << SDL_GetTicks() << "*** SDLFrame::nach Render" << endl;

			//     Point <int> p;
			//     RenderNamedScene(p);

			/*if (imgwriter && imgwriterflag==1) {
				imgwriter->Write(this);
				imgwriterflag=2;
			}*/

			/*  if (gl->screenshot==2) {
			  if (gl->iw)
			// MATTHIAS     gl->iw->Write();
			  gl->screenshot=1;
			 }
			 render = false;
			*/

	SDL_GL_SwapWindow(screen->window);
    	//cout << SDL_GetTicks() << "*** SDLFrame::nach Swapr" << endl;

		}

		if (Redraw()) render = true;


	  	//cout << SDL_GetTicks() << "*** SDLFrame::vor sleep" << endl;

		#ifdef WIN32
		Sleep(10);
#else

		usleep(3*1000);
#endif
  	//cout << SDL_GetTicks() << "*** SDLFrame::nach sleep" << endl;



		// if (SDL_SemValue(sem1))
		// SDL_SemPost(sem1);
	}


	SDL_RemoveTimer(my_timer_id);
	SDL_ShowCursor(SDL_ENABLE);
	return returnvalue;
}

/*
void SDLFrame::Info() {

 cout << "GL_RENDERER   = " << glGetString(GL_RENDERER) << endl;
 cout << "GL_VERSION    = " << glGetString(GL_VERSION) << endl;
 cout << "GL_VENDOR     = " << glGetString(GL_VENDOR) << endl;
 cout << "GL_EXTENSIONS = " << glGetString(GL_EXTENSIONS) << endl;

}
*/

void SDLFrame::accPerspective(GLdouble fovy, GLdouble aspect,
			      GLdouble nnear, GLdouble ffar, GLdouble pixdx,
			      GLdouble pixdy, GLdouble eyedx, GLdouble eyedy, GLdouble focus) {

	GLdouble fov2,left,right,bottom,top;
	fov2 = ((fovy*M_PI)/180.0)/2.0;
	top = nnear/(cos(fov2)/sin(fov2));
	bottom = -top;
	right = top*aspect;
	left = -right;
	accFrustum(left,right,bottom,top,nnear,ffar,pixdx,pixdy,eyedx,eyedy,focus);
}

void SDLFrame::accFrustum(GLdouble left, GLdouble right, GLdouble bottom,
			  GLdouble top, GLdouble nnear, GLdouble ffar, GLdouble pixdx,
			  GLdouble pixdy, GLdouble eyedx, GLdouble eyedy, GLdouble focus) {

	GLdouble xwsize, ywsize;
	GLdouble dx, dy;
	GLint viewport[4];

	glGetIntegerv(GL_VIEWPORT, viewport);

	xwsize = right - left;
	ywsize = top - bottom;
	dx = -(pixdx * xwsize / ((GLdouble) viewport[2]) + eyedx * nnear / focus);
	dy = -(pixdy * ywsize / ((GLdouble) viewport[3]) + eyedy * nnear / focus);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	{
		glFrustum(left + dx, right + dx, bottom + dy, top + dy, nnear, ffar);
	}

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(static_cast<GLfloat>(-eyedx), static_cast<GLfloat>(-eyedy), 0.0f);
}

bool SDLFrame::KeyboardFunc(SDL_Keysym&) {
	return false;
}
bool SDLFrame::SpecialKeys(const Uint8*) {
	return false;
}
bool SDLFrame::MouseButton(SDL_MouseButtonEvent&) {
	return false;
}
void SDLFrame::MouseWheel(const SDL_MouseWheelEvent&) {

}
bool SDLFrame::MouseMotion(SDL_MouseMotionEvent&) {
	return false;
}
bool SDLFrame::JoystickButtonEvent(SDL_JoyButtonEvent&) {
	return false;
}
bool SDLFrame::JoystickAxisEvent(SDL_JoyAxisEvent&) {
	return false;
}
bool SDLFrame::JoystickHatEvent(SDL_JoyHatEvent&) {
	return false;
}

void SDLFrame::TimerFunc() {}

void SDLFrame::Message(int, const void*) {}

void SDLFrame::Reshape(int w, int h) {
	screen->twwindow.size.x=w;
	screen->twwindow.size.y=h;
}

void SDLFrame::RenderScene() {}

int SDLFrame::Redraw() {
	return true;
}

}
