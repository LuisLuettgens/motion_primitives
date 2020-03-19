
#include "smoothmovement.h"


SmoothMovement::SmoothMovement() : value(0), only100(0), velocity(0), rot(false), gototime(0) {}


void SmoothMovement::Init(float minv, float maxv, bool r, int o100) {
	minval = minv;
	maxval = maxv;
	rot = r;
	only100 = o100;
	
}

void SmoothMovement::JoyInput(int j) {

	value = minval + (maxval-minval)*(((j/32767.f)+1.f)/2.f);

}


void SmoothMovement::Set(float v) {
	value = v;

	if (!rot) {
		if (value > maxval)
			value = maxval;

		if (value < minval)
			value = minval;
	}
}


bool SmoothMovement::Timer(double dt) {

	if (gototime>0) {

		value -= (value - gotoval)/gototime*dt;
		gototime -=dt;

		if (gototime<=0) {
			gototime=0;
			value = gotoval;
			return true;
		}

	}

	if (velocity) {

		value += dt * .005f * velocity * 1000.f;
		if (velocity > 0)
			velocity --;
		if (velocity < 0)
			velocity++;

		if (rot) {
			if (value > maxval)
				value -= maxval;

			if (value < minval)
				value += maxval;
		} else {
			if (value > maxval)
				value = maxval;

			if (value < minval)
				value = minval;
		}
	}

	if (only100 && velocity==0) {
		int v = (int)((value+50)/100);
		value = v*100.f;
	
	}

	return false;
}


void SmoothMovement::Accelerate(int step, int maxspeed) {
	/*if (velocity == 0)
		velocity = -step;
	velocity -= step;
	if (velocity < -maxspeed)
		velocity = -maxspeed;*/

	gototime = 0.;

	if (velocity == 0)
		velocity = step;
	velocity += step;
	if (maxspeed<0) {
		if (velocity < maxspeed)
			velocity = maxspeed;
	} else {
		if (velocity > maxspeed)
			velocity = maxspeed;
	}
}


void SmoothMovement::Goto(float v, double time) {

	if (rot) {
		while (v-value>180.f) value+=360.;
		while (v-value<-180.f) value-=360.;
	}

	gotoval = v;
	gototime = time;

}
