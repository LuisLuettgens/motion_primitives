#pragma once

#ifndef _WIN32
#include "TWGUIconfig.h"
#endif

#include "SDL2/SDL_thread.h"

#ifdef TRANSWORHP_GRAPHICS

namespace tw {

using threadf = void* (*) (int *active);

class TransWorhpProblem;
class TWfolder;

/** @ingroup gl
 *  @brief Starting jobs in different thread (e.g. optimizing).
 */
class SDLThread {
public:

	/** Constructor. */
	SDLThread();
	
	/** Destructor. */
	~SDLThread();
	
	/** Start job. */
	int Run(TWfolder *tw, TransWorhpProblem *ph);
	
	/** Lock data. */
	int Lock();
	
	/** Unlock data. */
	int Unlock();

	int threadlive;
	threadf threadfunction;
	
	int Active();
private:
	SDL_Thread *sdl_t;
	SDL_mutex *lock;
	SDL_sem *sem;
};

extern SDLThread *thethread0;

}

#endif
