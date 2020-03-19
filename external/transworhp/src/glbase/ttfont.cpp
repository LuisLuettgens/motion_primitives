#ifdef WIN32
#include "windows.h"
#endif

#include <GL/gl.h>
#include <string.h>
#include "ttfont.h"
#include "conversion.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include "textout.h"
#include "../imaging/imageformat.h"
#include "../base/twstatus.h"

using namespace std;

TTFontManager::TTFontManager() : W(32), H(32) {

	int error = FT_Init_FreeType( &library );

	if ( error ) {
		//cout << "WERWER WR" << error << endl;
	}

	FT_Int amajor, aminor, apatch;
	FT_Library_Version( library, &amajor, &aminor, &apatch );

	version = ToString(amajor) + "." + ToString(aminor) + "." + ToString(apatch);

	MyStatus("Version", string("Freetype ") + version, Status::NORMAL);
}

TTFontManager::~TTFontManager() {
	FT_Done_FreeType(  library );
}

TTFontData::~TTFontData() {

	glDeleteLists(fontOffset,255);



}

void TTFontManager::printStringCenter(int font, const char *s, float x, float y, float z) const {

	int l = StringLength(font,s);
	printString(font,s,x-l/2,y,z);
}

void TTFontManager::printString(int font, const char *s, float x, float y, float z) const {

	glPushMatrix();

	glEnable( GL_TEXTURE_2D );
	glTranslatef(x,y,z);

	glListBase(ttf[font].fontOffset);
	glCallLists(strlen(s), GL_UNSIGNED_BYTE, (GLubyte *) s);

	glDisable( GL_TEXTURE_2D );
	glPopMatrix();
}


void TTFontManager::printString(int font, const char *s, float x, float y, float z, int width) const {

	printString(font, s, x, y, z);

}


void TTFontManager::makeLetter(int font, int ch, const rectangle &tx, int dx) {

	//cout << "MAKELETTER " << ch << ":::" << endl;
	int c = ch;

	glNewList(ttf[font].fontOffset + ch, GL_COMPILE);

	if (c>128) {
		c-=128;
		ttf[font].tex2.Bind();
	} else {
		ttf[font].tex1.Bind();
	}
	glBegin(GL_QUADS);

	glTexCoord2f(tx.width/(float)W,(c+1)/128.f);
	glVertex3f((float)tx.width,0.f,0.f);

	glTexCoord2f((float)tx.width/(float)W,c/128.f);
	glVertex3f((float)tx.width,(float)tx.height,0.f);

	glTexCoord2f(0.0f,c/128.f);
	glVertex3f(0.f,(float)tx.height,0.f);

	glTexCoord2f(0.0f,(c+1)/128.f);
	glVertex3f(0,0,0.f);

	glEnd();

	glTranslatef(dx+0.f,0.f,.01f);

	glEndList();
}


void TTFontManager::Load(int font, const std::string &filename, int size) {

	FT_Face face;

	string st =	string("Loading font ") + filename + " @ " + ToString(size);

	int error = FT_New_Face( library, filename.c_str(), 0, &face );
	if ( error == FT_Err_Unknown_File_Format ) {
		//cout << " the font file " << filename
		//<< " could be opened and read, but it appears that its font format is unsupported" << endl;
		MyStatus("TrueType", st + " unsupported", Status::ERR);
		return;
	} else if ( error ) {
//		cout << " ... another error code means that the font file " << filename
//		<< " could not be opened or read, or simply that it is broken..." << endl;
		MyStatus("TrueType", st + " failed", Status::ERR);
		return;

	}
	MyStatus("TrueType", st + " ok", Status::WARN);

	error = FT_Set_Pixel_Sizes( face, /* handle to face object */
	                            0, /* pixel_width */
	                            size ); /* pixel_height */

	ttf[font].height = size;

	GLubyte *data = new GLubyte[W*H*256*4];
	for (int i=0;i<W*H*256*4;i++) {
		data[i] = 0;
	}

	for (int i=0;i<255;i++) {
		ttf[font].widths[i] = 0;
	}
	FT_GlyphSlot slot = face->glyph;

	int dp = 4;

	for (int charcode=0;charcode<255-1;charcode++) {

		unsigned char c = charcode;

		error = FT_Load_Char( face, c, FT_LOAD_RENDER );

		GLubyte *r = &data[charcode*W*H*dp];
		unsigned char *b = slot->bitmap.buffer;

		for (unsigned int y=0;y<slot->bitmap.rows;y++) {
			GLubyte *r2 = r + (y + 25-(slot->metrics.horiBearingY>>6) )*W*dp
			              + (slot->metrics.horiBearingX>>6) *dp;

			for (unsigned int x=0;x<slot->bitmap.width;x++) {
				*r2=255;
				r2++;
				*r2=255;
				r2++;
				*r2=255;
				r2++;
				*r2=*b;
				b++;
				r2++;
			}

		}

		/* increment pen position */
		ttf[font].widths[charcode] = (slot->advance.x >> 6 );
	}

	ttf[font].widths[(int)' '] = ttf[font].widths[(int)'n'];

	ttf[font].tex1.LoadTexture(data,0,W,H*128,dp);
	ttf[font].tex2.LoadTexture(&data[128*W*H*dp],0,W,H*128,dp);

	ttf[font].tex1.Bind();
	glEnable(GL_TEXTURE_2D);

	ttf[font].fontOffset = glGenLists(255);

	for (int charcode=31;charcode<255-1;charcode++) {

		unsigned char c = charcode;

		error = FT_Load_Char(face, c, FT_LOAD_RENDER);
		rectangle r(0,charcode*H,slot->bitmap.width+12,H);

		makeLetter(font,charcode,r,ttf[font].widths[charcode]);

	}
	glDisable(GL_TEXTURE_2D);

	delete []data;

FT_Done_Face(face);

}


int TTFontManager::StringLength(int font, const char *s) const {


	int ret=0;
	for (int i=0;s[i];i++) {

		int index = s[i];
		if (index<0)
			index+=256;

		ret += ttf[font].widths[index];

		//	cout << (char) index << index<< " "  << ttf[font].widths[index] << " " << ret << endl;
	}
	//cout << " " <<  ret << endl;
	return ret;
}

