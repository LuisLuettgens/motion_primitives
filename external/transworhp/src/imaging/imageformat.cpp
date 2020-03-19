#ifdef WIN32
#include "windows.h"
#endif
#include "imageformat.h"

#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

ImageFormat* ImageFormat::newFormat(const std::string &filename) {

	//std::cout << "newFormat(" << filename << ")" << std::endl;
	std::string suffix = filename.substr(filename.length()-3);
	//std::cout << suffix << std::endl;

	ImageFormat *ret=nullptr;
/*	if (suffix == "png") {
		ret = new PNGImageFormat(filename);
	} else if (suffix == "tif") {
		ret = new TIFFImageFormat(filename);
	} else if (suffix == "gif") {
		ret = new GIFImageFormat(filename);
	} else if (suffix == "bmp") {
		ret = new BMPImageFormat(filename);
	} else if (suffix == "jpg" || suffix == "jpeg" || suffix == "JPG" || suffix == "JPEG") {
		ret = new JPEGImageFormat(filename);
	}*/
	return ret;
}

ImageFormat::ImageFormat() : data(nullptr), type("") {}

ImageFormat::ImageFormat(const std::string &t) : data(nullptr), type(t) {}

ImageFormat::ImageFormat(const std::string &f_name, const std::string &t)
		: data(nullptr), filename(f_name), type(t) {}

ImageFormat::~ImageFormat() {
	delete [] data;
}

bool ImageFormat::Read() {
	return false;
}

bool ImageFormat::Write(GLubyte */*p*/) {
	return false;
}

color4 ImageFormat::Get(int x, int y) {

	color4 ret;
	if (data) {
		unsigned char *c = &data[(x+y*width)*depth];
		if (depth==4) {
			ret = color4(*c/256.f,*(c+1)/256.f,*(c+2)/256.f,*(c+3)/256.f);
		}
		if (depth==3) {
			ret = color4(*c/256.f,*(c+1)/256.f,*(c+2)/256.f);
		}
		if (depth==1) {
			ret = color4(1,1,1,*(c)/256.f);
		}
	}
	return ret;
}
