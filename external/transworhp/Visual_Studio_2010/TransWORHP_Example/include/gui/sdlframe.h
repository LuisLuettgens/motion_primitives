#ifndef sdlframe_h
#define sdlframe_h
#include <GL/gl.h>
#include "../base/defines.h"
 
#include <SDL2/SDL.h>
#include <string>
#include "xmlio.h"

#include "sdlscreen.h"
class ImageWriter;

#ifdef WIN32
#ifndef MINGW
#define DllExport __declspec( dllexport )
#pragma warning (disable : 4251)
#else
#define DllExport
#endif
#else 
#define DllExport
#endif

/** @ingroup gl
 *  @brief Create physical window under Linux and Windows with OpenGL context.
 */
class DllExport SDLFrame {
public:
	SDLFrame(SDLScreen* screen_);
	virtual ~SDLFrame();

	int Loop(int wait, int terminate);
	int Wait();
	void Timer();
	//void Info();

	void MyQuit(int status);

	int waittime;

SDLScreen *screen;


    virtual bool KeyboardFunc(SDL_Keysym &keysym) {
        return 0;
    }
    virtual bool SpecialKeys(const Uint8 *keys) {
        return 0;
    }
    virtual bool MouseButton(SDL_MouseButtonEvent &m) {
        return 0;
    }
    virtual bool MouseMotion(SDL_MouseMotionEvent &m) {
        return 0;
    }
    virtual bool JoystickButtonEvent(SDL_JoyButtonEvent &j) {
        return 0;
    }
    virtual bool JoystickAxisEvent(SDL_JoyAxisEvent &j) {
        return 0;
    }
    virtual bool JoystickHatEvent(SDL_JoyHatEvent &j) {
        return 0;
    }
    virtual void TimerFunc() {}


    virtual void Message(int id, const void *param=0) {}
    void Reshape(int w, int h) {screen->twwindow.size.x=w; screen->twwindow.size.y=h;}
    virtual void RenderScene() {};
    //virtual void RenderNamedScene()=0;


   
int returnvalue;


    void accPerspective(GLdouble fovy, GLdouble aspect,
                        GLdouble nnear, GLdouble ffar, GLdouble pixdx,
                        GLdouble pixdy, GLdouble eyedx, GLdouble eyedy, GLdouble focus);
    void accFrustum(GLdouble left, GLdouble right, GLdouble bottom,
                    GLdouble top, GLdouble nnear, GLdouble ffar, GLdouble pixdx,
                    GLdouble pixdy, GLdouble eyedx, GLdouble eyedy, GLdouble focus);
protected:



    ImageWriter *imgwriter;
    int imgwriterflag;
virtual int Redraw() {return true;}
  
private:

	bool render;
	int busy;

};

#endif

