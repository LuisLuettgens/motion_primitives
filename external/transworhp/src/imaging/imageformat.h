#pragma once

#include "../base/defines.h"
#include "../base/color4.h"
#include <GL/glew.h>
#include <GL/gl.h>
#include <string>

/** @defgroup imaging Load and Save images
 *  @brief *.png *.bmp *.tif *.gif
 */


/** @ingroup imaging
 *  @brief Interface for different image formats.
 *
 */
class ImageFormat {

public:
	/** Constructor. */
	ImageFormat();

	/** Constructor.
	 *
	 * @param t
	 */
	ImageFormat(const std::string &t);

	/** Constructor.
	 *
	 * @param f_name
	 * @param t
	 */
	ImageFormat(const std::string &f_name, const std::string &t);

	/** Destructor. */
	virtual ~ImageFormat();


	/** Load file.
	 *
	 * @return
	 */
	virtual bool Read();
	/** Save data.
	 *
	 * @param p
	 * @return
	 */
	virtual bool Write(GLubyte *p);


	/** Create instance of class according to file suffix.
	 *
	 * @param filename
	 * @return
	 */
	static ImageFormat *newFormat(const std::string &filename);

	color4 Get(int x, int y);

	int width, height, depth;
	unsigned char *data;

	std::string description;
	std::string filename;
	std::string type;
};


/** @ingroup imaging
 *  @brief Implementation for PNG. Uses libpng.
 *
 */
class PNGImageFormat : public ImageFormat {

public:
	PNGImageFormat();
	PNGImageFormat(const std::string &f_name);
	~PNGImageFormat();

	bool Read();
	bool Write(GLubyte *p);

};


/** @ingroup imaging
 *  @brief Implementation for TIFF. Uses libtif.
 *
 *  Only available if compiled with #define USE_TIFF
 *
 */
class TIFFImageFormat : public ImageFormat {

public:
	TIFFImageFormat();
	TIFFImageFormat(const std::string &f_name);
	~TIFFImageFormat();

	bool Read();
	bool Write(GLubyte *p);

};

/** @ingroup imaging
 *  @brief Implementation for GIF. Uses libungif.
 *
 *  Only available if compiled with #define USE_GIF
 *
 */
class GIFImageFormat : public ImageFormat {

public:
	GIFImageFormat();
	GIFImageFormat(const std::string &f_name);
	~GIFImageFormat();

	bool Read();
	bool Write(GLubyte *p);

};

/** @ingroup imaging
 *  @brief Implementation for JPEG. Uses libjpeg.
 *
 *  Only available if compiled with #define USE_JPEG
 *
 */

class JPEGImageFormat : public ImageFormat {

public:
	JPEGImageFormat();
	JPEGImageFormat(const std::string &f_name);
	~JPEGImageFormat();

	bool Read();
	bool Write(GLubyte *p);

private:
	int read ();

	void put_pixel_rows(void* cinfo, void* dinfo, int rows_supplied);
};

/** @ingroup imaging
 *  @brief Implementation for BMP.
 *
 */
class BMPImageFormat : public ImageFormat {

public:
	BMPImageFormat();
	BMPImageFormat(const std::string &f_name);
	~BMPImageFormat();

	bool Read();
	bool Write(GLubyte *p);

private:
	int bfType() {
		return get2(0);
	}
	int bfSize() {
		return get4(2);
	}
	int bfReserved1() {
		return get2(6);
	}
	int bfReserved2() {
		return get2(8);
	}
	int bfOffBits() {
		return get4(10);
	}
	int biSize() {
		return get4(14);
	}
	int biWidth() {
		return get4(18);
	}
	int biHeight() {
		return get4(22);
	}
	int biPlanes() {
		return get2(26);
	}
	int biBitCount() {
		return get2(28);
	}

	int get2(int a);
	int get4(int a);
	int set2(int a, int v);
	int set4(int a, int v);

	unsigned char header[54];
};
