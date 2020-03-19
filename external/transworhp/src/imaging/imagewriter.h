#pragma once

#include "../base/defines.h"
#include "xmlio.h"

#include <GL/glew.h>
#include <GL/gl.h>
#include <string>

namespace tw {

class SDLFrame;
class ToolWindow;

/** @ingroup imaging
 *  @brief Render image to file.
 *
 */
class ImageWriter {
public:
	/** Constructor. */
	ImageWriter(XMLNode *n);

	/** Destructor. */
	~ImageWriter();

	/** Render whole screen to file. */
	int Write(SDLFrame *gl);

	/** Render toolwindow to file. */
	int WriteTool(ToolWindow *tp,SDLFrame *gl);

private:
	int width,height;
	std::string filename;
	std::string description;

	int time;
	int index;

	GLubyte *image;
};

}
