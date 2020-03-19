//
// Author: Matthias Knauer <knauer@math.uni-bremen.de>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
#ifndef ttfont_h
#define ttfont_h
#include <GL/gl.h>

#include <string>
#include "xmlio.h"
#include "texture.h"
#include "font.h"
#include "viewport.h"

#include <ft2build.h>
#include FT_FREETYPE_H

struct TTFontData {
	TTFontData() {}

	~TTFontData();
	Texture tex1;
	Texture tex2;
	int widths[256];
	GLuint fontOffset;
	int height;

};

class TTFontManager {
public:
	TTFontManager();
~TTFontManager();

	void Load(int font, const std::string &filename, int size);

	void printString(int font, const char *s, float x, float y, float z) const;
	void printString(int font, const char *s, float x, float y, float z, int width) const;
	void printStringCenter(int font, const char *s, float x, float y, float z) const;
	int StringLength(int font, const char *s) const;

	TTFontData ttf[15];

	std::string version;
private:
	void makeLetter(int font, int ch, const rectangle &tx, int dx);
	FT_Library library;

	int W, H;

};




/*
class TTFont : public Font {
public:
	TTFont();
	virtual ~TTF();
 
	void makeLetter(char ch, const rectangle &tx);
	*/ /*
	protected:
		virtual void makeGLFont() {
			std::cout << "WER " << std::endl;
		}*/
/*
	Texture *tex;
};
 
*/




#endif
