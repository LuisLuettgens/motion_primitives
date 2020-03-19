#include "sdlthread.h"

#include <iostream>

#ifdef WIN32
#include "windows.h"
#endif

namespace tw {

SDLThread *thethread0 = nullptr;

int subthread(void *arg) {

	SDL_sem* f = (SDL_sem*)arg;

	while (true) {

		SDL_SemWait(f);

		if (thethread0) {
			if (thethread0->threadfunction) {
				thethread0->threadfunction(&thethread0->threadlive);
				//std::cout << "subthread" << std::endl;
				thethread0->threadlive = 2;
			} else {
				break;
			}
		} else {
			break;
		}
	}

	return 0;
}


SDLThread::SDLThread() : threadlive(1), threadfunction(nullptr), sdl_t(nullptr) {

	lock = SDL_CreateMutex();

	sem = SDL_CreateSemaphore(0);

	sdl_t = SDL_CreateThread(subthread, "WORHP Thread", (void*) sem);
}


int SDLThread::Active() {
	
	int a = SDL_SemValue(sem);
	if (!a) {
		return 0;
	}
	else {
		std::cout << "SDLThread still active: " << a << std::endl;
		return 1;
	}
}

SDLThread::~SDLThread() {

//cout << "Delete SDLThread" << endl;
	threadlive = 0;

	int a = SDL_SemValue(sem);
	if (!a) {
		threadfunction = nullptr;
		SDL_SemPost(sem);
	} else {
		std::cout << "SDLThread still active: " << a << std::endl;
	}

	threadfunction = nullptr;

	if (sdl_t) {
		SDL_WaitThread(sdl_t, nullptr);
	}

	Lock();
	Unlock();

	SDL_DestroyMutex(lock);
	SDL_DestroySemaphore(sem);
}


TransWorhpProblem *globalTransWorhp = nullptr;
TWfolder *globaltwfolder = nullptr;

void* calcTransWorhp(int *);

int SDLThread::Run(TWfolder *tw, TransWorhpProblem *ph) {

	globalTransWorhp = ph;
	globaltwfolder = tw;

	if (SDL_SemValue(sem) > 0) {
		return 0;
	} else {
		threadfunction = calcTransWorhp;
		//cout << "MATT RUN" << endl;
		SDL_SemPost(sem);
	}
	return 0;
}


int SDLThread::Lock() {
	return SDL_mutexP(lock);
}


int SDLThread::Unlock() {
	return SDL_mutexV(lock);
}

}
