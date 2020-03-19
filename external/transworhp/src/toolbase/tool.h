#pragma once

#include "../glbase/texture.h"
#include "../glbase/font.h"

#include "../base/point.h"
#include "../base/color4.h"

#include <SDL2/SDL.h>

namespace tw {

class SDLFrame;
/**
 * Implement this function:
 * Send message.
 * @param id Message ID.
 */
void Message(int id);

/**
 * Implement this function:
 * Send message with argument.
 * @param id Message ID.
 * @param param Argument.
 */
void Message(int id, const void *param);

/**
 * Implement this function:
 * Add status information.
 * @param text Status information.
 */
///void Status(const std::string &text);


typedef void (*CallbackFkt) ();


/** @ingroup toolbase
 *  @brief Tool.
 *
 * Static variables and common routines for OpenGL Tool System.
 *
 */
class Tool {
public:

	/**
	 * Initialize OpenGL Tool System.
	 * @param xml XML Node with TOOL tag.
	 */
	static void InitStatic(XMLNode *xml);

	/**
	 * Draw Rectangle
	 *
	 * @param p1 Lower left point.
	 * @param p2 Upper right point.
	 * @param z  Depth.
	 * @param mode Rectangle style. OR combination of these flags:
	 *        - 1: Draw inner region bright.
	 *        - 2: Draw inner region dark.
	 *        - 4: Draw rectangle inset.
	 *        - 8: Draw rectangle outset.
	 *        - 16: Draw rectangle with dark border.
	  - 256: Center.
	 */
	static void drawRect(const Point<int> &p1,const Point<int> &p2,float z, int mode);

	/**
	 * Convert char in SDLKey.
	 * @param c Character.
	 * @return SDLKey.
	 */
	static SDL_Keycode toKey(char c);

	/** Texture for graphical elements. */
	static Texture *texture;

	/** Standard font. */
	static FFont *font;

	/** Set of used colors. */
	static color4 colors[10];

	/** Active Mode, depending on last selected window. */
	static int activeMode;

	static std::string filename;
	static std::string filename_prefix;

	enum DrawMode_e {FILL=1, FILL2=2, INSET=4, OUTSET=8, SOLID=16, FILLDARK=32, FILLHI=64, CENTER=256};

	static SDLFrame *frame;
};

template <class T>
void clearpointer(std::vector<T *> &vlist) {
	typename std::vector<T *>::iterator it = vlist.begin();
	for (;it!=vlist.end();it++) {
		delete *it;
	}
	vlist.clear();
}

}
