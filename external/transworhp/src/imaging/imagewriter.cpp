#include "imagewriter.h"

#include "conversion.h"

#include "../imaging/imageformat.h"

#include "../toolbase/toolwindow.h"

#include "../gui/sdlthread.h"
#include "../gui/sdlframe.h"

#include "../core/twstatus.h"

#include <iostream>

using namespace std;

namespace tw {

ImageWriter::ImageWriter(XMLNode *n) {

	filename = n->GetText();

	/*std::string s;
	s = n->GetAttribute("width");
	if (s!="")
	 width = ToInt(s);

	s = n->GetAttribute("height");
	if (s!="")
	 height = ToInt(s);
	*/

	//compression = COMPRESSION_NONE;

	//description= n->GetAttribute("description");

	time =0;
	index=0;

	width = 2000;
	height = 2000;
	image = new GLubyte[width * height * 3];

}


ImageWriter::~ImageWriter() {

	delete []image;
}


int ImageWriter::WriteTool(ToolWindow *tp, SDLFrame *gl) {

	int drawbuffer[1];

	glGetIntegerv(GL_DRAW_BUFFER, &drawbuffer[0]);

	glDrawBuffer(GL_BACK);

	int x=0, y=0;

	width = tp->width;
	height = tp->height;

	int glwidth = gl->screen->Width();
	int glheight = gl->screen->Height();

	if (glwidth<width || glheight<height)
		return 1;

	width = (width/4)*4;
	height = (height/4)*4;

	gl->Reshape(width,height);
	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);

	//glClear(GL_COLOR_BUFFER_BIT || GL_DEPTH_BUFFER_BIT); -> BITWISE OR
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	if (tp)  {
		thethread0->Lock();
		tp->RenderWindow();
		//gl->BackgroundColor();
		thethread0->Unlock();
	}


	glGetIntegerv(GL_VIEWPORT, viewport);

	glReadBuffer(GL_BACK);
	glReadPixels(x, y, width, height, GL_RGB, GL_UNSIGNED_BYTE, image);

	char buf[100];
	sprintf(buf,"%s_%s.png",Tool::filename_prefix.c_str(),tp->GetFilename().c_str());

	ImageFormat *img=ImageFormat::newFormat(std::string(buf));

	if (img) {
		img->width=width;
		img->height=height;
		img->description = description;

		if (img->Write(image)) {
			//cout << "Writing " << buf << endl;

			stringstream st;
			st << buf;
			MyStatus("Writing", st.str(), Status::NORMAL);
			// index++;
		}

		delete img;
	}

	glDrawBuffer(drawbuffer[0]);
	glReadBuffer(drawbuffer[0]);
	gl->Reshape(glwidth,glheight);
	gl->RenderScene();

	return 0;
}


int ImageWriter::Write(SDLFrame *gl) {
	/*
	 int numBuffers[10];
	 glGetIntegerv(GL_AUX_BUFFERS, &numBuffers[0]);
	*/
	//cout << "AUX=" <<numBuffers[0]<< endl;
	int drawbuffer[1];

	glGetIntegerv(GL_DRAW_BUFFER, &drawbuffer[0]);

	//glDrawBuffer(GL_AUX0);
	glDrawBuffer(GL_BACK);

	int x=0, y=0;

	width = gl->screen->Width();
	height = gl->screen->Height();
	//std::cout << width << " " << height << std::endl;

	width = (width/4)*4;
	height = (height/4)*4;
	//std::cout << width << " " << height << std::endl;

	int ww = gl->screen->Width();
	int hh = gl->screen->Height();
	gl->Reshape(width,height);

	//GLuint BlurTexture = EmptyTexture();
	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	//cout << viewport[0] << " " << viewport[1] << " " << viewport[2] << " " << viewport[3] << endl;

	//glClearColor(0,0.0,0.3,0);
	//glClear(GL_COLOR_BUFFER_BIT || GL_DEPTH_BUFFER_BIT); -> BITWISE OR
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

//((glBahn*)gl)->hidemenu=0;
	gl->RenderScene();
//((glBahn*)gl)->hidemenu=1;

	glGetIntegerv(GL_VIEWPORT, viewport);

	//glBindTexture(GL_TEXTURE_2D,BlurTexture);

	//glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_LUMINANCE, 0, 0, 1024, 1024, 0);


	//glReadBuffer(GL_AUX0);
	glReadBuffer(GL_BACK);
	// glPixelStorei(GL_PACK_ALIGNMENT, 1);

	glReadPixels(x, y, width, height, GL_RGB, GL_UNSIGNED_BYTE, image);

	char buf[100];
	sprintf(buf,filename.c_str(),index);

	ImageFormat *img=ImageFormat::newFormat(std::string(buf));

	if (img) {
		img->width=width;
		img->height=height;
		img->description = description;
		//cout << img->width << endl;

		if (img->Write(image)) {
			//cout << "Writing " << buf << endl;
			index++;

			stringstream st;
			st << buf;
			MyStatus("Image", st.str(), Status::NORMAL);
		}

		delete img;

	}

	glDrawBuffer(drawbuffer[0]);
	glReadBuffer(drawbuffer[0]);
	gl->Reshape(ww,hh);
	//gl->RenderScene(true);

	return 0;
}

}
