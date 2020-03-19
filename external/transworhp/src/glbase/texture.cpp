/***************************************************************************
*   Copyright (C) 2004 by Matthias Knauer                                 *
*   knauer@math.uni-bremen.de                                             *
*                                                                         *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or     *
*   (at your option) any later version.                                   *
*                                                                         *
*   This program is distributed in the hope that it will be useful,       *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
*   GNU General Public License for more details.                          *
*                                                                         *
*   You should have received a copy of the GNU General Public License     *
*   along with this program; if not, write to the                         *
*   Free Software Foundation, Inc.,                                       *
*   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
***************************************************************************/
#ifdef WIN32
#include "windows.h"
#endif
#include "texture.h"
#include <iostream>
#include <fstream>
#include "../imaging/imageformat.h"
#include "../base/exception.h"


#define GL_CLAMP_TO_EDGE                  0x812F


Texture::Texture() : texture(0), data(0), width(0), height(0), dim(0) {dataw=0;datah=0; }


Texture::~Texture() {
	glDeleteTextures( 1, &texture );

	if (data)
		delete []data;
}

int Texture::GetData(int x, int y, int d) const {
	return data[(x+width*y)*depth+d];
}

void Texture::Bind() const {

	if (texture && dim==2)
		glBindTexture(GL_TEXTURE_2D, texture);

	if (texture && dim==1)
		glBindTexture(GL_TEXTURE_1D, texture);

}

int Texture::GetWidth() const {
	return width;
}
int Texture::GetHeight() const {
	return height;
}


// load a 256x256 RGB .RAW file as a texture
void Texture::LoadTextureRAW( const char * filename, int wrap, int w, int h, int dp) {

	dim = 2;
	width = w;
	height = h;
	depth = dp;
	if (h==1)
		dim=1;

	std::ifstream file( filename,std::ios::binary );

	//	std::cout << "Loading " << filename << std::endl;
	// open texture data
	if ( file.bad() )
		return;

	// allocate buffer
	GLubyte *dx = new GLubyte[width * height * depth];

	// read texture data
	char a;
	// Soeren: Matthias bist du sicher, dass die obere Schleife sinnvoll ist?
	// Entweder 54 <= width * height * depth => dx[i] darf geschrieben werden,
	// aber dann überschreibt die Schleife direkt danach die Werte sowieso wieder.
	// Wenn 54 > width * height* depth => Speicherzugriffsfehler. Mir unklar was hier passieren soll
	for (int i=0;i<54;i++) {
		file.get(a);
		dx[i] = a;
	}
	for (int i=0;i<width * height * depth;i++) {
		file.get(a);
		dx[i] = a;
	}


	// Soeren: Matthias auch hier sagt die Codeanalyse, dass die Zuweisung: dx[i*depth] = dx[i*depth+2];
	// Probleme verursachen könnte.
	for (int i=0;i<width*height;i++) {
		char tmp=dx[i*depth];
		dx[i*depth] = dx[i*depth+2];
		dx[i*depth+2] = tmp;
	}

	LoadTexture(dx,wrap,w,h,dp);

}


// load a 256x256 RGB .RAW file as a texture
void Texture::LoadTexture( GLubyte *d, int wrap, int w, int h, int dp ) {
	int error;

	data = d;
	dim = 2;
	width = w;
	height = h;
	depth = dp;
	if (h==1)
		dim=1;

	// allocate a texture name
	glGenTextures( 1, &texture );
error = glGetError();

	int GLTEX = GL_TEXTURE_2D;
	int GLRGB = GL_RGB;

	if (dim==1)
		GLTEX = GL_TEXTURE_1D;
	if (dp==1) {
		GLRGB = GL_RED;
		glPixelTransferf(GL_RED_SCALE,1./256.);
		glPixelTransferf(GL_GREEN_SCALE,1./256.);
		glPixelTransferf(GL_BLUE_SCALE,1./256.);
		glPixelTransferf(GL_ALPHA_SCALE,1./256.);
		glPixelTransferf(GL_RED_BIAS,1./256.);
		glPixelTransferf(GL_GREEN_BIAS,1.);
		glPixelTransferf(GL_BLUE_BIAS,0.);
		glPixelTransferf(GL_ALPHA_BIAS,0.);
	}
	if (dp==4)
		GLRGB = GL_RGBA;

	// select our current texture
	glBindTexture( GLTEX, texture );
	error = glGetError();

	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	error = glGetError();



	// if wrap is true, the texture wraps over at the edges (repeat)
	//       ... false, the texture ends at the edges (clamp)
	glTexParameterf( GLTEX, GL_TEXTURE_WRAP_S,
		wrap ? GL_REPEAT : GL_CLAMP );
	error = glGetError();
	glTexParameterf( GLTEX, GL_TEXTURE_WRAP_T,
		wrap ? GL_REPEAT : GL_CLAMP );
	error = glGetError();

	glTexParameteri(GLTEX, GL_TEXTURE_MAG_FILTER,
		GL_LINEAR);
	error = glGetError();
	glTexParameteri(GLTEX, GL_TEXTURE_MIN_FILTER,
		GL_LINEAR);
	error = glGetError();

	// select modulate to mix texture with color for shading
	glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE );
	error = glGetError();
	//glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);

	float col[4]= {
		1,1,1,1
	};
	glTexParameterfv( GLTEX, GL_TEXTURE_BORDER_COLOR,col);
	error = glGetError();


	if (dim==2)
		glTexImage2D(GL_TEXTURE_2D, 0, depth, width,
		height,	0, GLRGB, GL_UNSIGNED_BYTE,
		data);
	else

		glTexImage1D(GL_TEXTURE_1D, 0, depth, width,
		0, GLRGB, GL_UNSIGNED_BYTE,
		data);

	error = glGetError();
	if (error) std::cout << "glTexImage2D unsupported parameter: " << error << std::endl;
	//int A =  GL_MAX_TEXTURE_SIZE;


	/*
	// when texture area is small, bilinear filter the closest mipmap
	glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
	GL_LINEAR_MIPMAP_NEAREST );
	// when texture area is large, bilinear filter the first mipmap
	glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
	*/
	// build our texture mipmaps
	//// MATTHIAS
	///	   gluBuild2DMipmaps( GL_TEXTURE_2D, 3, width, height, GL_RGB, GL_UNSIGNED_BYTE, data );
	//

	/*   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexImage2D(GL_TEXTURE_2D, 0, 3, width, height, 0, GL_RGB , GL_UNSIGNED_BYTE, data);
	*/

	//delete  []data ;

	if (dp==1) {
		glPixelTransferf(GL_RED_SCALE,1.);
		glPixelTransferf(GL_GREEN_SCALE,1.);
		glPixelTransferf(GL_BLUE_SCALE,1.);
		glPixelTransferf(GL_ALPHA_SCALE,1.);
		glPixelTransferf(GL_RED_BIAS,0.);
		glPixelTransferf(GL_GREEN_BIAS,0.);
		glPixelTransferf(GL_BLUE_BIAS,0.);
		glPixelTransferf(GL_ALPHA_BIAS,0.);
	}
	
	data = 0;
}



void Texture::LoadTexturePNG(const std::string &filename, int wrap, bool auto_alpha ) {




	ImageFormat *img = 0;
	bool ret=0;

	if (!data) {
		width=0;
		height=0;
		depth = 4;
		dim=2;

		img = ImageFormat::newFormat(filename);
		if (!img)
			return;

		ret = img->Read();
		//img.GetData(data);
	}

	if (ret) {

		format = GL_RGBA;
		if (img->depth==3)
			format = GL_RGB;
		//std::cout << "D" << img->depth << std::endl;
		width=img->width;
		height=img->height;

		//	std::cout << img->width << std::endl;
		//	std::cout << img->height << std::endl;

		if (height<=0 || width<=0)
			return;

		/*if (height==1 || height==2 || height==4 || height==8 || height==16 || height==32
			|| height==64 || height==128 || height==256 || height==512
			|| height==1024 || height==2048 || height==4096) {}
		else {
			Exception e;
			e << "Image " << filename << ": height=" << height << " is not a power of 2";
			throw e;
		}
		if (width==1 || width==2 || width==4 || width==8 || width==16 || width==32
			|| width==64 || width==128 || width==256 || width==512
			|| width==1024 || width==2048 || width==4096) {}
		else {
			Exception e;
			e << "Image " << filename << ": width=" << width << " is not a power of 2";
			throw e;
		}*/


		if (img->type=="PNG") {
			depth = img->depth;

			data = new GLubyte[width*height*depth];
			for (int i=0;i<width*height*depth;i++) {
				data[i] = img->data[i];
			}
		} else if (img->type=="JPG") {
			depth = img->depth;

			//std::cout << width << " " << height << " " << depth << std::endl;
			data = new GLubyte[width*height*depth];
			for (int i=0;i<width*height*depth;i++) {
				data[i] = img->data[i];
				//std::cout << data[i] << ":" ;
			}
		} else if (img->type=="BMP") {


			if (auto_alpha == false) {
				if (img->depth!=24) {
					Exception e;
					e << "Image " << filename << ": depth=" << img->depth << " is not 24";
					throw e;

				}
				depth = 3;  //img.depth/8; //3
				format = GL_RGB;

				data = new GLubyte[width*height*depth];
				for (int i=0;i<width*height*depth;i++) {
					data[i] = img->data[i];
				}

			} else {
				if (img->depth!=24) {
					Exception e;
					e << "Image " << filename << ": depth=" << img->depth << " is not 24";
					throw e;
				}
				depth = 4;

				data = new GLubyte[width*height*depth];
				for (int i=0;i<width*height;i++) {
					data[i*depth] = img->data[i*3];
					data[i*depth+1] = img->data[i*3+1];
					data[i*depth+2] = img->data[i*3+2];
					data[i*depth+3] = 255-(img->data[i*3]+img->data[i*3+1]+img->data[i*3+2])/3;
				}

			}

		}

	}

	//std::cout << "LOAD " << filename << " " << data[0]  << " " << width << " " << height<< std::endl;
	if (data) {
		// allocate a texture name
		glGenTextures( 1, &texture );

		//	std::cout << "OUT" << texture << std::endl;
		// select our current texture
		glBindTexture( GL_TEXTURE_2D, texture );

		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

		// if wrap is true, the texture wraps over at the edges (repeat)
		//       ... false, the texture ends at the edges (clamp)
		glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S,
			wrap ? GL_REPEAT : GL_CLAMP_TO_EDGE );
		glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T,
			wrap ? GL_REPEAT : GL_CLAMP_TO_EDGE );

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,
			GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
			GL_LINEAR);

		// select modulate to mix texture with color for shading
		glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE );
		//glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);

		float col[4]={1,1,1,1};
		glTexParameterfv( GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR,col);

		glTexImage2D(GL_TEXTURE_2D, 0, depth, width,
			height, 0, format, GL_UNSIGNED_BYTE, data);


	}

	delete img;

#ifndef WIN32
	//	delete []data;
	//	data = 0 ;
#endif

}



void Texture::GrabScreen(int dx, int dy, int w, int h) {

	if (dataw*datah < w*h) {
		delete []data;
		data=0;
	}

	if (data==0) {

		format = GL_RGBA;
		depth = 4;
		width = w;
		height = h;
		dim = 2;
dataw=w;
datah=h;

		data = new GLubyte[width*height*depth];


		if (data==0) {
			fprintf(stderr,"Failed to allocate memory for the texture\n");
			return;
		}

		glReadPixels(dx,dy,width,height,format,GL_UNSIGNED_BYTE,data);

		glGenTextures(1,&texture);
		glBindTexture(GL_TEXTURE_2D,texture);
		glTexParameterf(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,GL_CLAMP);
		glTexParameterf(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,GL_CLAMP);
		glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
		glTexImage2D(GL_TEXTURE_2D,0,4, width, height, 0,format,GL_UNSIGNED_BYTE,data);
		glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_MODULATE);

	}
	else {
		glReadPixels(dx,dy,width,height,format,GL_UNSIGNED_BYTE,data);
		glBindTexture(GL_TEXTURE_2D,texture);
		glTexImage2D(GL_TEXTURE_2D,0,4, width, height,  0,format,GL_UNSIGNED_BYTE,data);
	}

}