#ifdef WIN32
#include "windows.h"
#endif
#include "imageformat.h"
#include <GL/gl.h>
#include <iostream>
#include <fstream>

#include "png.h"
#include <sstream>
#include "../base/status.h"
#include "conversion.h"
#include <string.h>

using namespace std;

PNGImageFormat::PNGImageFormat() : ImageFormat("PNG") {}

PNGImageFormat::PNGImageFormat(const std::string &f_name) : ImageFormat(f_name,"PNG") {}

PNGImageFormat::~PNGImageFormat() {}


void user_read_data(png_structp png_ptr, png_bytep data, png_size_t length) {

	/*voidp read_io_ptr = png_get_io_ptr(png_ptr);
	std::ifstream *file = (std::ifstream*) read_io_ptr;
	file->read((char*)data,length);*/
}

void user_write_data(png_structp png_ptr, png_bytep data, png_size_t length) {

	/*voidp write_io_ptr = png_get_io_ptr(png_ptr);
	std::ofstream *file = (std::ofstream*) write_io_ptr;
	file->write((char*)data, length);*/
}

void user_flush_data(png_structp png_ptr) {

	/*voidp write_io_ptr = png_get_io_ptr(png_ptr);
	std::ofstream *file = (std::ofstream*) write_io_ptr;
	file->flush();*/
}



bool PNGImageFormat::Read() {

	png_image image;

	/* Only the image structure version number needs to be set. */
	memset(&image, 0, sizeof image);
	image.version = PNG_IMAGE_VERSION;

	if (png_image_begin_read_from_file(&image, filename.c_str())) {
		png_bytep buffer;

		/* Change this to try different formats!  If you set a colormap format
		* then you must also supply a colormap below.
		*/
		image.format = PNG_FORMAT_RGBA;

		buffer = (png_bytep) malloc(PNG_IMAGE_SIZE(image));

		if (buffer != NULL) {
			if (png_image_finish_read(&image, NULL/*background*/, buffer,
				0/*row_stride*/, NULL/*colormap for PNG_FORMAT_FLAG_COLORMAP */))
			{


				width = image.width;
				height = image.height;
				depth = 4; // image.png_get_channels(png_ptr, info_ptr);


//png_bytepp row_pointers;
//row_pointers = png_get_rows(png_ptr, info_ptr);



data = new unsigned char[width*height*depth];
for (int i=0; i<height; i++) {
for (int j=0; j<width*depth; j++) {
data[j+width*depth*i] = buffer[j+width*depth*i]; // row_pointers[i][j];
}
}



				free(buffer);
			}

			else
			{
				MyStatus("PNG", "pngtopng: read " + filename + ": " + image.message,
				         Status::NORMAL);

				/* This is the only place where a 'free' is required; libpng does
				* the cleanup on error and success, but in this case we couldn't
				* complete the read because of running out of memory.
				*/
				png_image_free(&image);
				return false;
			}
		}

		else {
			MyStatus("PNG", string("pngtopng: out of memory: ") + to_string(PNG_IMAGE_SIZE(image)) +
			                    "bytes",
			         Status::NORMAL);
			return false;
		}
	}

	else {
		/* Failed to read the first argument: */

		MyStatus("PNG", "pngtopng: " + filename + ": " + image.message, Status::NORMAL);
		return false;
	}

  return true;
}











/* Set error handling if you are using the setjmp/longjmp method (this is
* the normal method of doing things with libpng).  REQUIRED unless you
* set up your own error handlers in the png_create_read_struct() earlier.
*/
/*	if (setjmp(png_jmpbuf(png_ptr))) {
/* Free all of the memory associated with the png_ptr and info_ptr */
/*		png_destroy_read_struct(&png_ptr, &info_ptr, png_infopp_NULL);
return false;
}

/* Set up the input control if you are using standard C streams */
//png_init_io(png_ptr, fp);
/*	png_set_read_fn(png_ptr, &file, user_read_data);




/* If we have already read some of the signature */
/*	png_set_sig_bytes(png_ptr, sig_read);

png_read_png(png_ptr, info_ptr, 0, 0);


/* At this point you have read the entire image */


/*	width = png_get_image_width(png_ptr, info_ptr);
height = png_get_image_height(png_ptr, info_ptr);

//int bit_depth = png_get_bit_depth(png_ptr, info_ptr);
//int color_type = png_get_color_type(png_ptr, info_ptr);
int channels = png_get_channels(png_ptr, info_ptr);
//int rowbytes = png_get_rowbytes(png_ptr, info_ptr);

depth = channels;

//std::cout << "R" << rowbytes << std::endl;
//std::cout << width << " " << height << " " << bit_depth
//<< " " << color_type << " " << channels << std::endl;

png_bytepp row_pointers;
row_pointers = png_get_rows(png_ptr, info_ptr);



data = new unsigned char[width*height*channels];
for (int i=0; i<height; i++) {
for (int j=0; j<width*channels; j++) {
data[j+width*channels*i] = row_pointers[i][j];
}
}

/* clean up after the read, and free any memory allocated - REQUIRED */
/*	png_destroy_read_struct(&png_ptr, &info_ptr, png_infopp_NULL);


/* that's it */



/*
void PNGImageFormat::GetData(GLubyte *d) {

// geht nicht....
d = new GLubyte[width*height*4];
for (int i=0;i<width*height*4;i++) {
std::cout << (int) data[i] << " ";
d[i] = data[i];

}


//png_bytepp row_pointers;
//row_pointers = png_get_rows(png_ptr, info_ptr);


//d = new GLubyte[width*height*channels];
//for (int i=0; i<height; i++) {
//for (int j=0; j<width*channels; j++) {
// d[j+width*channels*i] = row_pointers[i][j];
//}
//}

}*/


bool PNGImageFormat::Write(GLubyte *p) {

	std::ofstream file(filename.c_str(), std::ios::binary);
	// std::clog << "Loading " << filename << std::endl;

	if (file.bad()) {
		std::cerr << "Failed writing " << filename << std::endl;
		return false;
	}


	png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, 0, 0, 0);
	if (!png_ptr) {
		return false;
	}


	png_infop info_ptr = png_create_info_struct(png_ptr);
	if (!info_ptr) {
		png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
		return false;
	}


	if (setjmp(png_jmpbuf(png_ptr))) {
		png_destroy_write_struct(&png_ptr, &info_ptr);
		return false;
	}


	png_set_write_fn(png_ptr, &file, user_write_data, user_flush_data);



	int bit_depth = 8;
	int color_type = PNG_COLOR_TYPE_RGB;
	int interlace_type = PNG_INTERLACE_NONE;
	int compression_type = PNG_COMPRESSION_TYPE_DEFAULT;
	int filter_method = PNG_FILTER_TYPE_DEFAULT;

	png_set_IHDR(png_ptr, info_ptr, width, height,
		bit_depth, color_type, interlace_type,
		compression_type, filter_method);


	float gamma = .5;
	png_set_gAMA(png_ptr, info_ptr, gamma);


	/* Optionally write comments into the image */
	png_text text_ptr[3];
	char desc[100];
	strcpy(desc,description.c_str());
	text_ptr[0].key = (char*)"Title";
	text_ptr[0].text = desc;
	text_ptr[0].compression = PNG_TEXT_COMPRESSION_NONE;
	text_ptr[1].key = (char*)"Author";
	text_ptr[1].text = (char*)"Matthias Knauer";
	text_ptr[1].compression = PNG_TEXT_COMPRESSION_NONE;
	text_ptr[2].key = (char*)"Description";
	text_ptr[2].text = desc;
	text_ptr[2].compression = PNG_TEXT_COMPRESSION_NONE;
#ifdef PNG_iTXt_SUPPORTED

	text_ptr[0].lang = NULL;
	text_ptr[1].lang = NULL;
	text_ptr[2].lang = NULL;
#endif

	png_set_text(png_ptr, info_ptr, text_ptr, 3);


	png_bytep* row_pointers = new png_bytep[height];
	for (int i=0;i<height;i++) {
		row_pointers[height-1-i] = &p[i*width*3];
	}


	/// &row???
	png_set_rows(png_ptr, info_ptr, row_pointers);


	png_write_png(png_ptr, info_ptr, 0, 0);


	/* clean up after the read, and free any memory allocated - REQUIRED */
	png_destroy_write_struct(&png_ptr, &info_ptr);

	delete []row_pointers;
	return true;
}



