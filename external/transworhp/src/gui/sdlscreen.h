#pragma once

#include "../base/point.h"

#include "xmlio.h"

#include <SDL2/SDL.h>


#include <GL/glew.h>
#include <GL/gl.h>

namespace tw {

class SDLFrame;
class SDLThread;

struct TWwindow {
	
	TWwindow();
	
	/** XML lesen.
	 * @param xmlmain XML-Knote
	 * @param countParams zum Parameter zeahlen: Gesamtzahl aller Parameter
	 * @param setParams #gesetzter Parameter
	 */
	void ParseXML(XMLNode *xml, int *countParams=nullptr, int *setParams=nullptr);
	
	/** Drawable OpenGL Area */
	Point<int> size;
	
	/** Physical screen resolution */
	Point<int> resolution;
	
	/** Displacement of upper left corner */
	Point<int> reference;
	
	bool fullscreen;
	bool maximized;
	
	int multisamplebuffers;
	int multisamplesamples;
	
	/** 0: normal, 1: quiet */
	int windowmode;
};



/** @ingroup gl
 *  @brief Create physical screen under Linux and Windows.
 */
class SDLScreen {
public:
	SDLScreen(TWwindow *twscreen);
	virtual ~SDLScreen();
	
	// creates shader
	void initGL();

	void ToggleFullScreen(SDLFrame *parent);
	void SetIcon();
	
	int Width() const;
	int Height() const;

	SDLThread *thethread;
	SDL_Surface *image;
	
	SDL_Window *window;
	SDL_GLContext glcontext;

	TWwindow twwindow;
	
	GLuint simpleProgram;
	GLuint textureProgram;
	GLuint sparsePlotProgram;
};

}
