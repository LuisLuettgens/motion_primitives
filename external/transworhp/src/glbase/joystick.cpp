#include "joystick.h"

#include <iostream>
#include "conversion.h"
#include <sstream>
#include "../base/twstatus.h"

using namespace std;



Joystick::Joystick() : joy(0), numhats(0) {

	// Initialize the joystick subsystem
	SDL_InitSubSystem(SDL_INIT_JOYSTICK);

	index = -1;
	int num_joy = SDL_NumJoysticks();

{
	stringstream ss;
	ss << num_joy << (num_joy==1?" joystick":" joysticks") << " found";
	MyStatus("Joystick",ss.str(),num_joy?2:4);
}


	for (int i=0;i<num_joy;i++) {
		string s(SDL_JoystickNameForIndex(i));

		stringstream ss;
		ss << "Name: " << s << " ";

		if (s == "Saitek Saitek Pro Flight Yoke" || s == "Saitek Pro Flight Yoke") {
			name = s;
			index = i;
			ss << "*";
		}

    /*    if (strncmp("Generic",s.c_str(),7)==0) {
			name = s;
			index = i;
			ss << "*";
		}
*/
		MyStatus("Joystick", ss.str(), Status::NORMAL);
	}

	// Check for joystick
	if(index!=-1) {
		// Open joystick
		joy=SDL_JoystickOpen(index);

		if (joy) {
			//	printf("Opened Joystick 0\n");
			//	printf("Name: %s\n", SDL_JoystickName(index));

			MyStatus("Joystick", "Number of Axes: " + ToString(SDL_JoystickNumAxes(joy)),
			         Status::NORMAL);
			MyStatus("Joystick", "Number of Buttons: " + ToString(SDL_JoystickNumButtons(joy)),
			         Status::NORMAL);
			MyStatus("Joystick", "Number of Balls: " + ToString(SDL_JoystickNumBalls(joy)),
			         Status::NORMAL);
			MyStatus("Joystick", "Number of Hats: " + ToString(SDL_JoystickNumHats(joy)),
			         Status::NORMAL);

			/*
			            printf("Number of Axes: %d\n", SDL_JoystickNumAxes(joy));
			            printf("Number of Buttons: %d\n", SDL_JoystickNumButtons(joy));
			            printf("Number of Balls: %d\n", SDL_JoystickNumBalls(joy));
			            printf("Number of Hats: %d\n", SDL_JoystickNumHats(joy));*/
			numhats = SDL_JoystickNumHats(joy);

			MyStatus("Joystick", string("Opened Joystick ") + name + " #" + ToString(index),
			         Status::WARN);
		} else
			MyStatus("Joystick", string("Couldn't open Joystick #" + ToString(index)), Status::ERR);

	} else {
		;//MyStatus("Joystick", "No joystick found!");
	}

	SDL_JoystickEventState(SDL_ENABLE);

}


Joystick::~Joystick() {

	// Close if opened
	if(SDL_JoystickOpen(index))
		SDL_JoystickClose(joy);

}

int Joystick::GetAxis(int axis) {

	return SDL_JoystickGetAxis(joy, axis);

}

bool Joystick::GetButton(int button) {

	return SDL_JoystickGetButton(joy, button);

}

int Joystick::GetHat(int hat) {

	if (numhats) {
		return SDL_JoystickGetHat(joy, hat);
	} else {
		// fake a hat for linux

		int ret = 0;

		int xm = GetAxis(5);
		int ym = GetAxis(6);

		if (xm<0)
			ret |= SDL_HAT_LEFT;
		if (xm>0)
			ret |= SDL_HAT_RIGHT;
		if (ym<0)
			ret |= SDL_HAT_UP;
		if (ym>0)
			ret |= SDL_HAT_DOWN;

		return ret;

	}
}

/*
void Joystick::ButtonEvent(SDL_JoyButtonEvent &e) {

	cout << (int)e.type << " : " << (int)e.which << " : "
	<< (int)e.button << " : " << (int)e.state << endl;

}
*/
