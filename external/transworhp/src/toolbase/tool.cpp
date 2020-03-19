#include "tool.h"

#ifdef WIN32
#include "windows.h"
#endif

namespace tw {

Texture *Tool::texture = nullptr;
FFont *Tool::font = nullptr;
int Tool::activeMode = 0;
std::string Tool::filename = "polyearth.xml";
std::string Tool::filename_prefix = "polyearth";
SDLFrame* Tool::frame = nullptr;

color4 Tool::colors[10] = {
                              color4(.95f, .7f, .7f, .9f),
                              color4(.95f, .93f, .93f, .9f),
                              color4(0.f, 0.f, 0.f, 1.f)
                          };


void Tool::InitStatic(XMLNode *xml) {

	// Titelleiste
	std::string s = xml->GetAttribute("color_back");
	if (s!="") {
		colors[0] = s;
	}

	// Hauptfeld
	s = xml->GetAttribute("color_extra");
	if (s!="") {
		colors[6] = s;
	}

	// Schrift im Titel
	s = xml->GetAttribute("color_fps");
	if (s!="") {
		colors[5] = s;
	}

	// Kontrastfarbe
	/*s = xml->GetAttribute("color4");
	if (s!="") {
		colors[8] = s;
	}
	else {
		const float *a = colors[0].GetData();
		const float *b = colors[6].GetData();
		float ab = .75;
		colors[8] = color4(ab*a[0]+(1-ab)*b[0],ab*a[1]+(1-ab)*b[1],ab*a[2]+(1-ab)*b[2],ab*a[3]+(1-ab)*b[3]);
	}*/

	colors[2] = colors[0].Bright(.75f);
	colors[1] = colors[0].Dim(.5f);
	colors[3] = colors[0].Bright(.9f);
	colors[4] = colors[0].Dim(.8f).Alpha(1.0f);
	colors[7] = colors[0].Dim(.15f);


}



// 1 : inneres hell
// 2: inneres dunkel
// 4 : rein
// 8 : erhaben
// 16 : dunkel

void Tool::drawRect(const Point<int> &p1, const Point<int> &p2, float z, int mode) {

	if (mode & (1+2+FILLDARK+FILLHI)) {
		if (mode&2) {
			glColor4fv(colors[0].GetData());
		} else if (mode&FILLDARK) {
			glColor4fv(colors[7].GetData());
		} else if (mode&FILLHI) {
			glColor4fv(colors[7].GetData());
		} else {
			glColor4fv(colors[2].GetData());
		}

		glBegin(GL_QUADS);

		glVertex3f(static_cast<GLfloat>(p1.x), static_cast<GLfloat>(-p1.y), z);
		glVertex3f(static_cast<GLfloat>(p1.x), static_cast<GLfloat>(-p2.y), z);

		glVertex3f(static_cast<GLfloat>(p2.x), static_cast<GLfloat>(-p2.y), z);
		glVertex3f(static_cast<GLfloat>(p2.x), static_cast<GLfloat>(-p1.y), z);
		glEnd();
	}

	if (mode & (4+8+16)) {
		glColor4fv(colors[mode&8?3:1].GetData());
		glBegin(GL_LINE_STRIP);
		glVertex3f(static_cast<GLfloat>(p2.x+1), static_cast<GLfloat>(-p1.y), z+.1f);
		glVertex3f(static_cast<GLfloat>(p1.x), static_cast<GLfloat>(-p1.y), z+.1f);
		glVertex3f(static_cast<GLfloat>(p1.x), static_cast<GLfloat>(-p2.y), z+.1f);
		glEnd();

		glColor4fv(colors[mode&4?3:1].GetData());
		glBegin(GL_LINE_STRIP);
		glVertex3f(static_cast<GLfloat>(p1.x), static_cast<GLfloat>(-p2.y), z+.1f);
		glVertex3f(static_cast<GLfloat>(p2.x), static_cast<GLfloat>(-p2.y), z+.1f);
		glVertex3f(static_cast<GLfloat>(p2.x), static_cast<GLfloat>(-p1.y), z+.1f);
		glEnd();
	}
}


SDL_Keycode Tool::toKey(char c) {

	SDL_Keycode h = SDLK_CLEAR;
	if (c>='a' && c<='z') {
		h = (SDL_Keycode)(c + SDLK_a - 'a');
	} else if (c>='A' && c<='Z') {
		h = (SDL_Keycode)(c + SDLK_a - 'A');
	} else if (c>='0' && c<='9') {
		h = (SDL_Keycode)(c + SDLK_0 - '0');
	} else {
		switch (c) {
		case '.':
			h=SDLK_PERIOD;
			break;
		case ',':
			h=SDLK_COMMA;
			break;
		case '-':
			h=SDLK_MINUS;
			break;
		case '#':
			h=SDLK_HASH;
			break;
		}
	}

	return h;
}

}
