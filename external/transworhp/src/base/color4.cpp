#include "color4.h"

#include "conversion.h"

#include <cmath>
#include <cstdio>
#include <iostream>

struct colorname {
	char name[20];
	int c[3];
};

colorname colorlist[]  = {
                             {"Black",{0,0,0}},
                             {"Red",{255,0,0}},
                             {"Green",{0,255,0}},
                             {"Blue",{0,0,255}},
                             {"Yellow",{255,255,0}},
                             {"Magenta",{255,0,255}},
                             {"Cyan",{0,255,255}},
                             {"White",{255,255,255}},

                             {"Purple",{77,17,202}},
                             {"PoisonPurple",{77,5,102}},
                             {"Violet",{100,100,200}},
                             {"ForestGreen",{77,178,102}},
                             {"FreshGreen",{77,178,12}},
                             {"LimeGreen",{177,217,2}},
                             {"Orange",{227,178,12}},
                             {"Sand",{177,178,102}},
                             {"Ochre",{177,178,2}},
                             {"Steel",{77,78,102}},
                             {"SteelBlue",{77,117,202}},
                             {"Salmon",{177,117,102}},
                             {"Coral",{177,78,102}},
                             {"Coral2",{200,100,100}},

                             {"SkyBlue",{240,240,255}},
                             {"Bronze",{155,121,88}},
                             {"",{255,255,255}} // Ende
                         };


color4::color4() : name("White") {
	c[0] = 1.;
	c[1] = 1.;
	c[2] = 1.;
	c[3] = 1.;
}

color4::color4(float r, float g, float b, float a) {
	c[0] = r;
	c[1] = g;
	c[2] = b;
	c[3] = a;
}

color4::color4(const std::string &a) {

	if (a[0]>='0' && a[0]<='9') {

		std::vector<int> aa = ToIntArray(a);
		c[0] = aa[0]/255.f;
		c[1] = aa[1]/255.f;
		c[2] = aa[2]/255.f;
		c[3] = 1.f;

		if (aa.size()>=4)
			c[3] = aa[3]/255.f;
	} else {
		std::vector<std::string> aa = ToStringArray(a);
		int i = 0;
		while (true) { //colorlist[i].name!="") {
			if (aa[0] == colorlist[i].name || colorlist[i].name[0]==0) {
				c[0] = colorlist[i].c[0]/255.f;
				c[1] = colorlist[i].c[1]/255.f;
				c[2] = colorlist[i].c[2]/255.f;

				name = colorlist[i].name;
				break;
			}
			i++;
		}

		if (aa.size()>=2)
			c[3] = ToInt(aa[1])/255.f;
		else
			c[3] = 1.f;

	}
}


void color4::SetGrey(float val) {
	c[0] = val;
	c[1] = val;
	c[2] = val;
	c[3] = 1.;
}


void color4::SetHue(float hue) {

	if (hue<1.f/5.f) {
		c[0] = 1.f;
		c[1] = hue*5.f;
		c[2] = 0.f;
	} else if (hue<2.f/5.f) {
		c[0] = 2.f - hue*5.f;
		c[1] = 1.f;
		c[2] = 0.f;
	} else if (hue<3.f/5.f) {
		c[0] = 0.f;
		c[1] = 1.f;
		c[2] = hue*5.f-2.f;
	} else if (hue<4.f/5.f) {
		c[0] = 0.f;
		c[1] = 4.f-hue*5.f;
		c[2] = 1.f;
	} else if (hue<1.) {
		c[0] = hue*5.f-4.f;
		c[1] = 0.f;
		c[2] = 1.f;
	} else {
		c[0] = 1.f;
		c[1] = 0.f;
		c[2] = 1.f;
	}
	c[3] = 1.f;
}


std::string ToString(const color4 &v) {
	if (v.name!="") {
		if (v.c[3] == 1.f) {
			return v.name;
		}
		return v.name + " " + ToString((int)(v.c[3] * 255.f));
	} else {
		char buf[30];
		sprintf(buf, "%d,%d,%d,%d",(int)(v.c[0]*255.f),(int)(v.c[1]*255.f),(int)(v.c[2]*255.f),
		        (int)(v.c[3]*255.f));
		return std::string(buf);
	}
}


color4 color4::Dim(float a) const {
	return color4((1.f-a)*c[0]+0*a, (1.f-a)*c[1]+0*a, (1.f-a)*c[2]+0*a,c[3]);
}

color4 color4::Bright(float a) const {
	return color4((1.f-a)*c[0]+1*a, (1.f-a)*c[1]+1*a, (1.f-a)*c[2]+1*a,c[3]);
}

color4 color4::Alpha(float a) const {
	return color4(c[0], c[1], c[2], c[3]*a);
}

const float* color4::GetData() const {
	return &c[0];
}

float color4::alpha() const {
	return c[3];
}
float color4::red() const {
	return c[0];
}
float color4::green() const {
	return c[1];
}
float color4::blue() const {
	return c[2];
}

bool color4::operator==(const color4 &other) const {
	if (c[0] != other.c[0]) { return false; }
	if (c[1] != other.c[1]) { return false; }
	if (c[2] != other.c[2]) { return false; }
	if (c[3] != other.c[3]) { return false; }
	return true;
}


std::ostream& operator<<(std::ostream &os, const color4 &c) {
	os << "(" << c.c[0] << "," << c.c[1] << "," << c.c[2] << ")";
	return os;
}


#ifdef TRANSWORHP_GRAPHICS
#include <GL/glew.h>
#include <GL/gl.h>

void color4::operator()() const {
	glColor4fv(c);
}

#endif
