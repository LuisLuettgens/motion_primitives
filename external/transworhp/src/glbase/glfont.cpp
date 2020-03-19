#include "glfont.h"

#include "tex_font.h"

static
#include "xmlfont2.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>

GLFont::GLFont() {

	XMLParser parser;
	XMLNode *node = parser.ParseString(xmlfile);

	if (node) {
		XMLNode *n = parser.GetDocument();
		makeGLFont(n,font2);
	} else {
		DebugStream d(standard_textoutputtype);
		d << beginfile;
		parser.GetError(d);
		d << *parser.GetDocument();
		d << endfile;
		std::cout << d.GetString();
	}
}


GLFont::GLFont(const std::string &filename) {
	XMLParser parser;
	XMLNode *node = parser.Parse(filename);

//	std::cout << "Loading font " << filename << std::endl;

	if (node) {
		XMLNode *n = parser.GetDocument();
		makeGLFont(n,font2);
	} else {
		DebugStream d(standard_textoutputtype);
		d << beginfile;
		parser.GetError(d);
		d << *parser.GetDocument();
		d << endfile;
		std::cout << d.GetString();
	}
}

GLFont::~GLFont() {

	glDeleteLists(fontOffset,255);
}


void GLFont::printString(const char *s, float x, float y, float z) const {

	glPushMatrix();
	tex.Bind();
	glEnable( GL_TEXTURE_2D );
	glTranslatef(x,y,z);

	glListBase(fontOffset);
	glCallLists(strlen(s), GL_UNSIGNED_BYTE, (GLubyte *) s);

	glDisable( GL_TEXTURE_2D );
	glPopMatrix();
}

void GLFont::printString(const char *s, float x, float y, float z, int /*width*/) const {

	//if (StringLength(s)<width) {
		printString(s,x,y,z);
	/*} else {
		char buf[100];
		strncpy(buf,s,100);
		buf[strlen(buf)-1]=0;
		while (StringLength(buf)>=width || buf == 0) {
			buf[strlen(buf)-1]=0;
		}
		printString(buf,x,y,z);

	}*/
}


void GLFont::makeGLFont(XMLNode *node, GLubyte *f) {

	XMLNode *n = node->GetFirstChild(std::string("LETTER"));
	while (n) {

		letter a(n);
		letters.push_back(a);

		n = node->GetNextChild(std::string("LETTER"));
	}

	tex.LoadTexture(f,0,256,256,4);

	tex.Bind();
	glEnable(GL_TEXTURE_2D);

	fontOffset = glGenLists (255);

	for (unsigned char ch=1; ch<255; ch++) {
		const letter *s = Key(ch);

		if (s) {

			float a = 1.f/16.f;
			float a1 = (ch%16)*a;
			float a2 = (ch/16)*a;
			float dax = 1.f/256.f*(s->width-1);
			float day = 1.f/256.f*12;
			int dy = -1;
			glNewList(fontOffset + ch, GL_COMPILE);

			if (s->character!=' ') {
				glBegin(GL_QUADS);
				glTexCoord2f(a1+dax,a2+day);
				glVertex3f(s->width-1,0+dy,0.f);

				glTexCoord2f(a1+dax,a2);
				glVertex3f(s->width-1,12+dy,0.f);

				glTexCoord2f(a1,a2);
				glVertex3f(0.f,12+dy,0.f);

				glTexCoord2f(a1,a2+day);
				glVertex3f(0.f,0+dy,0.f);

				glEnd();
			}
			glTranslatef(s->width,s->height,0.f);

			//glBitmap(8, 10, 0.0, 0.0, s->width, s->height, s->raster);
			//std::cout << ch << s->character << " " << s->width << std::endl;
			glEndList();

		}
	}


	glDisable(GL_TEXTURE_2D);
}
