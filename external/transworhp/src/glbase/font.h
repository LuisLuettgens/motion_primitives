#pragma once

#include "xmlio.h"
#include <GL/glew.h>
#include <GL/gl.h>
#include <string>

/** @ingroup glbase
 *  @brief Interface for OpenGL Fonts.
 *
 */
class FFont {
public:
	/** Konstruktor. */
	FFont() = default;
	/** Destruktor. */
	virtual ~FFont() = default;

	/** Ausgabe der Zeichenkette s an (x/y/z). */
	virtual void printString(const char *s, float x, float y, float z) const = 0;
	/** Ausgabe der Zeichenkette s an (x/y/z). */
	virtual void printString(const char *s, float x, float y, float z, int width) const = 0;

	/** Bestimmung der Breite der Zeichenkette. */
	virtual int StringLength(const char *s) const;


	static const unsigned char psi, phi, lambda,alpha,beta,gamma,delta,epsilon,Omega,omega;
	static const unsigned char UP, DOWN, STEP, BACK, DOT, UPUP, DOWNDOWN,CHECK,UNCHECK,GORIGHT;

	struct letter {
		letter(const std::string &text);
		letter(XMLNode *node);
		unsigned char character;
		int width;
		int height;
		GLubyte raster[10];
	};

	const letter* Key(unsigned char s) const;

protected:
	GLuint fontOffset;
	std::vector<letter> letters;
};
