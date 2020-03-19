#pragma once

#include <GL/glew.h>
#include <GL/gl.h>

#include "texture.h"
#include "font.h"

#include "xmlio.h"

#include <string>

/** @ingroup glbase
 *  @brief Smooth font, texture based.
 *
 */
class GLFont : public FFont {
public:
	GLFont();
	GLFont(const std::string &filename);

	~GLFont();

	void printString(const char *s, float x, float y, float z) const;
	void printString(const char *s, float x, float y, float z, int width) const;

private:
	void makeGLFont(XMLNode *node, GLubyte *f);
	Texture tex;
};
