#include "light.h"

#include "xmlio.h"

#include "conversion.h"

using namespace std;

Light::Light() : light(0) {}



void Light::Init(GLenum l, XMLNode *xml) {
	light = l;
	if (light == GL_LIGHT0) {
		ambientLight = color4(0.3f,0.3f,0.3f,1.f);
		diffuseLight = color4(0.6f,0.6f,0.6f,1.0f);
		specularLight = color4(1.0f,1.0f,1.0f,1.0f);
		spotCutoff = 0.;
		spotExponent = 0.;
		SetPos(1, 2., 5.);

	} else if (light == GL_LIGHT1) {
		ambientLight = color4(0.f,0.f,0.f,1.f);
		diffuseLight = color4(6.f,1.f,1.f,1.0f);
		specularLight = color4(1.0f,1.0f,1.0f,1.0f);
		spotCutoff = 40.;
		spotExponent = 90.;
	}
	else if (light == GL_LIGHT2) {
		ambientLight = color4(0.f,0.f,0.f,0.f);
		diffuseLight = color4(6.f,6.f,6.f,1.0f);
		specularLight = color4(1.0f,1.0f,1.0f,1.0f);
		spotCutoff = 20.;
		spotExponent = 28.;
	}

	if (xml) {
		string p = xml->GetAttribute("pos");
		if (p!="") {
			vector<double> pp = ToDoubleArray(p);
			SetPos(pp[0],pp[1],pp[2]);
			pos0[0] = pp[0];
			pos0[1] = pp[1];
			pos0[2] = pp[2];
			pos0[3] = 0;

		}

 p = xml->GetAttribute("direction");
		if (p!="") {
			vector<double> pp = ToDoubleArray(p);
			//SetPos(pp[0],pp[1],pp[2]);
			direction[0] = pp[0];
			direction[1] = pp[1];
			direction[2] = pp[2];
			direction[3] = pp[3];
		}

		p = xml->GetAttribute("ambient");
		if (p!="") {
			ambientLight = color4(p);
		}

		p = xml->GetAttribute("diffuse");
		if (p!="") {
			diffuseLight = color4(p);
		}

		p = xml->GetAttribute("specular");
		if (p!="") {
			specularLight =color4(p);
		}

		p = xml->GetAttribute("spotcutoff");
		if (p!="") {
			spotCutoff = (float)ToDouble(p);
		}
		p = xml->GetAttribute("spotexponent");
		if (p!="") {
			spotExponent = ToDouble(p);
		}
	}

	// Setup and enable light 0
	Set();

	if (light == GL_LIGHT0)
		Enable();

}

void Light::Enable() {
	if (light)
		glEnable(light);
}

void Light::Disable() {
	if (light)
		glDisable(light);
}

void Light::Position() {
	if (light)
		glLightfv(light, GL_POSITION, position);
}
void Light::Direction() {
	if (light)
		glLightfv(light, GL_SPOT_DIRECTION, direction);
}

void Light::Set() {
	if (light) {
		glLightfv(light, GL_AMBIENT, ambientLight.GetData());
		glLightfv(light, GL_DIFFUSE, diffuseLight.GetData());
		glLightfv(light, GL_SPECULAR, specularLight.GetData());

		if (spotCutoff) {
			glLightf(light, GL_SPOT_CUTOFF, spotCutoff);
			glLightf(light, GL_SPOT_EXPONENT, spotExponent);
		}
	}
}
