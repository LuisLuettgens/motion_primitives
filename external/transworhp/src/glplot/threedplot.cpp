#include "threedplot.h"

#include "xopt_eps.h"

#include "../glbase/globject.h"

#include "../gui/sdlframe.h"

#include "../toolbase/tool.h"

#include "../xmlio/conversion.h"

#include <iostream>
#include <iomanip>


using namespace std;

namespace tw {

void Camera::Init(XMLNode *n) {

	if (n) {
		string s( n->GetAttribute("pos") );
		if (!s.empty()) {
			vector<float> zz = ToFloatArray(s);
			x = zz[0];
			y = zz[1];
			z = zz[2];
		}
		s = n->GetAttribute("rotx");
		if (!s.empty()) rotx= ToFloat(s);

		s = n->GetAttribute("roty");
		if (!s.empty()) roty = ToFloat(s);

		s = n->GetAttribute("focus");
		if (!s.empty()) focus = ToFloat(s);

	}
}

/*
void ThreeDPlot::displayObjects(double *x, double t) {

	func3d(obj[0],x,t);

}*/

extern float jj[3][8][2];


ThreeDPlot::ThreeDPlot(XMLNode *xml, plot3d f) : BasePlot(0) {

	func3d = f;
	campos[0].Init(-1000.f,1000.f,false);
	campos[1].Init(-1000.f,1000.f,false);
	campos[2].Init(-1000.f,1000.f,false);
	campos[3].Init(0.f,360.f,true);
	campos[4].Init(-90.f,90.f,false);
	campos[5].Init(0.f,1000.f,false);

	string path("./");

	if (xml) {
		XMLNode *n = xml->GetFirstChild("PATH");
		if (n) {
			path = n->GetText();
		}

		n = xml->GetFirstChild("MODEL");
		while (n) {
			auto g = new glObject();

			g->Init(n);
			obj.push_back(g);
			//g->InitObject(path);
			g->InitSplitObject(path);

			n = xml->GetNextChild("MODEL");
		}

		n = xml->GetFirstChild("CAMERA");
		int i=0;
		while (n) {
			cam[i].Init(n);
			i++;
			n = xml->GetNextChild("CAMERA");
		}

		campos[0].Set(cam[0].x);
		campos[1].Set(cam[0].y);
		campos[2].Set(cam[0].z);
		campos[3].Set(cam[0].rotx);
		campos[4].Set(cam[0].roty);
		campos[5].Set(cam[0].focus);
	}
}



ThreeDPlot::~ThreeDPlot() {

	for (auto &glObj : obj) {
		delete glObj;
	}
	
	obj.clear();
}


void ThreeDPlot::Draw(DataStorage &/*ds*/, DataStorage &/*dstop*/) {

	Point<int> pos = GetPos();

	glPushMatrix();
	glTranslatef(static_cast<float>(pos.x),static_cast<float>(pos.y),0.f);

	//if (IsIcon())
	// drawIconFrame();
	//else
	drawBack();

	if (IsIcon()) {
		//drawIconFrame();
	} else {
		drawPlotName();
	}
	glPopMatrix();
}


void ThreeDPlot::drawBack() const {

	int rborder = 0;
	int lborder = 0;
	int bborder = 0;
	int tborder = 0;

	int width = GetWidth();
	int height = GetHeight();

	glBegin(GL_QUADS);
	float a = GetHue()*.2f+.3f; // [.3 ; .5]

	//int z1 = (int)((width-rborder-lborder)*.2);
	//int z2 = (int)((height-tborder-bborder)*.2);

	//int z = min(z1,z2);
	//z = min(z1,10);

	glColor4f(1.f,1.f,1.f,a);
	glVertex3f((float)lborder,(float)tborder,0.f);
	glVertex3f((float)(width-rborder),(float)tborder,0.f);
	glColor4f(1.f,1.f,1.f,a+.4f);
	glVertex3f(width-rborder,height-bborder,0.f);
	glVertex3f(lborder,height-bborder,0.f);
	glColor4f(1.f,1.f,1.f,a);
	glVertex3f(lborder,tborder,0.f);
	glEnd();

	glBegin(GL_LINE_STRIP);
	glColor4f(0.f,0.f,0.f,.5f);
	glVertex3f(lborder,tborder,0.1f);
	glVertex3f(width-rborder,tborder,0.1f);
	glVertex3f(width-rborder,height-bborder,0.1f);
	glVertex3f(lborder,height-bborder,0.1f);
	glVertex3f(lborder,tborder,0.1f);
	glEnd();

	glBegin(GL_LINE_STRIP);
	glColor4f(0.f,0.f,0.f,.1f);
	glVertex3f(0.f,0.f,0.1f);
	glVertex3f(0.f,height,0.1f);
	glVertex3f(width,height,0.1f);
	glVertex3f(width,0.f,0.1f);
	glVertex3f(0.f,0.f,0.1f);
	glEnd();
}





void ThreeDPlot::DrawMore(SDLFrame *viewer, double *x, double t) {

	GLint viewport0[4];
	glGetIntegerv(GL_VIEWPORT, viewport0);

	GLfloat ambientLight[] = { 0.3f, 0.3f, 0.3f, 1.0f };
	GLfloat diffuseLight[] = { 0.7f, 0.7f, 0.7f, 1.0f };

	glEnable(GL_DEPTH_TEST);	// Hidden surface removal
	//glFrontFace(GL_CCW);		// Counter clock-wise polygons face out
	glEnable(GL_CULL_FACE);		// Do not calculate inside of jet

	// Enable lighting
	glEnable(GL_LIGHTING);

	// Setup and enable light 0
	glLightfv(GL_LIGHT0,GL_AMBIENT,ambientLight);
	glLightfv(GL_LIGHT0,GL_DIFFUSE,diffuseLight);
	glEnable(GL_LIGHT0);
	/*
		// Enable color tracking
		glEnable(GL_COLOR_MATERIAL);

		// Set Material properties to follow glColor values
		glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
		*/
	//GLfloat nRange = 150.0f;
	GLfloat lightPos[] = { -50.f, 50.0f, 100.0f, 1.0f };


	// Set Viewport to window dimensions

	int rborder = 0;
	int lborder = 1;
	int bborder = 0;
	int tborder = 0;

	glViewport(viewport0[2]/2+ geom.cur.pos.x+lborder,viewport0[3]/2+ geom.cur.pos.y+rborder,
	           geom.cur.width-rborder-lborder, geom.cur.height-tborder-bborder);

	// Reset coordinate system
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glLightfv(GL_LIGHT0,GL_POSITION,lightPos);




	//GLfloat fovy=60;

	// Clear the window with current clearing color
	//glClear(GL_COLOR_BUFFER_BIT);

	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);

	/*glEnable(GL_ALPHA_TEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);*/
	int acsize1=1;

	if (acsize1==1) {
		viewer->accPerspective(campos[5].value,(GLdouble)viewport[2]/(GLdouble)viewport[3],
		                       .1,200.,jj[0][0][0],jj[0][0][1],0.,0.,1.);

		glPushMatrix();

		glTranslatef(campos[0].value,-campos[1].value,campos[2].value);
		glRotatef(campos[4].value, 1.0f, 0.0f, 0.0f);
		glRotatef(campos[3].value, 0.0f, 1.0f, 0.0f);

		glRotatef(-90.f, 1.0f, 0.0f, 0.0f);

		//glTranslatef(0,0,-70);
		//glRotatef( 10, 0.0f, 1.0f, 0.0f);
		//glRotatef(-80, 1.0f, 0.0f, 0.0f);

		if (obj.size()) {
			func3d(obj[0],x,t);
		} else {
			func3d(nullptr,x,t);
		}
		glPopMatrix();


	} else {
		int jitter;

		//int jindex = 2;
		//if (acsize1==2) jindex=0;
		//if (acsize1==4) jindex=1;

		for (jitter=0;jitter<acsize;jitter++) {
			glClear(GL_DEPTH_BUFFER_BIT);

			//	viewer->accPerspective(fovy,(GLdouble)viewport[2]/(GLdouble)viewport[3],
			//	               1.,4.,jj[jindex][jitter][0],
			//		       jj[jindex][jitter][1],0.,0.,1.);




			if (jitter==0)
				glAccum(GL_LOAD,1./acsize1);
			else
				glAccum(GL_ACCUM,1./acsize1);
		}

		glAccum(GL_RETURN,1.);
		glFlush();

		//	glClear(GL_DEPTH_BUFFER_BIT);

		//	glMatrixMode( GL_MODELVIEW );
		//	glLoadIdentity();


	}


	/*
	glDisable(GL_DEPTH_TEST);	// Hidden surface removal
		glFrontFace(GL_CW);		// Counter clock-wise polygons face out
		glDisable(GL_CULL_FACE);

		* /
		glDisable(GL_ALPHA_TEST);
		glEnable(GL_BLEND);
		/ *
		//glDisable(GL_STENCIL_TEST);


		glEnable(GL_COLOR_MATERIAL);
		glDisable(GL_LIGHTING);
		*/

	//glDisable(GL_TEXTURE);
	/*	glDisable(GL_TEXTURE_2D);
		glEnable(GL_COLOR_MATERIAL);*/

	glDisable(GL_LIGHTING);

}




bool ThreeDPlot::SpecialKeys(const Uint8 *keys) {


	int step = 3;
	int maxspeed = 25;
	/*



	else if (keys[SDLK_RSHIFT]) {
	if (keys[SDLK_UP])FlightTime+=1;
	if (keys[SDLK_DOWN])FlightTime-=1;
	//cout << "FlightTime "<< FlightTime << endl;
	}
	*/

	/*if (keys[SDLK_LSHIFT]) {
	    if (keys[SDLK_UP])
	        lightpos[2].Accelerate(step,maxspeed);
	    if (keys[SDLK_DOWN])
	        lightpos[2].Accelerate(-step,-maxspeed);
	    if (keys[SDLK_RIGHT])
	        lightpos[0].Accelerate(-step,maxspeed);
	    if (keys[SDLK_LEFT])
	        lightpos[0].Accelerate(step,-maxspeed);
	    if (keys[SDLK_PAGEUP])
	        lightpos[1].Accelerate(step,maxspeed);
	    if (keys[SDLK_PAGEDOWN])
	        lightpos[1].Accelerate(-step,-maxspeed);

	}
	else */
	if (keys[SDL_SCANCODE_RSHIFT]) {

		if (keys[SDL_SCANCODE_DOWN])
			campos[2].Accelerate(step,maxspeed);
		if (keys[SDL_SCANCODE_UP])
			campos[2].Accelerate(-step,-maxspeed);
		if (keys[SDL_SCANCODE_LEFT])
			campos[0].Accelerate(-step,maxspeed);
		if (keys[SDL_SCANCODE_RIGHT])
			campos[0].Accelerate(step,-maxspeed);
		if (keys[SDL_SCANCODE_PAGEDOWN])
			campos[1].Accelerate(step,maxspeed);
		if (keys[SDL_SCANCODE_PAGEUP])
			campos[1].Accelerate(-step,-maxspeed);
	}
	/*}
	else
	* /

	else if (keys[SDLK_RCTRL]) {

	    if (keys[SDLK_LEFT])
	        eyedist.Accelerate(-step,maxspeed);
	    if (keys[SDLK_RIGHT])
	        eyedist.Accelerate(step,-maxspeed);


	}
	*/
	else {
		if (keys[SDL_SCANCODE_LEFT])
			campos[3].Accelerate(-step,maxspeed);
		if (keys[SDL_SCANCODE_RIGHT])
			campos[3].Accelerate(step,-maxspeed);
		if (keys[SDL_SCANCODE_UP])
			campos[4].Accelerate(-step,maxspeed);
		if (keys[SDL_SCANCODE_DOWN])
			campos[4].Accelerate(step,-maxspeed);
		if (keys[SDL_SCANCODE_PAGEUP])
			campos[5].Accelerate(-step,maxspeed);
		if (keys[SDL_SCANCODE_PAGEDOWN])
			campos[5].Accelerate(step,-maxspeed);


	}/*


	if (keys[SDLK_w]) {
	    double ang = -campos[3].value*M_PI/180;
	    campos[2].Accelerate(cos(ang)*step,cos(ang)*maxspeed);
	    campos[0].Accelerate(sin(ang)*step,sin(ang)*maxspeed);
	}
	if (keys[SDLK_a]) {
	    double ang = campos[3].value*M_PI/180;
	    campos[2].Accelerate(sin(ang)*step,sin(ang)*maxspeed);
	    campos[0].Accelerate(cos(ang)*step,cos(ang)*maxspeed);
	}
	if (keys[SDLK_s]) {
	    double ang = -(campos[3].value+180)*M_PI/180;
	    campos[2].Accelerate(cos(ang)*step,cos(ang)*maxspeed);
	    campos[0].Accelerate(sin(ang)*step,sin(ang)*maxspeed);
	}
	if (keys[SDLK_d]) {
	    float ang = campos[3].value*M_PI/180;
	    campos[2].Accelerate(-sin(ang)*step, -sin(ang)*maxspeed);
	    campos[0].Accelerate(-cos(ang)*step, -cos(ang)*maxspeed);
	}*/






	/*
	else if (keys[SDLK_LCTRL]) {
	if (keys[SDLK_UP])
	cp[1].Accelerate(step,maxspeed);
	if (keys[SDLK_DOWN])
	cp[1].Accelerate(-step,-maxspeed);
	if (keys[SDLK_LEFT])
	cp[0].Accelerate(-step,maxspeed);
	if (keys[SDLK_RIGHT])
	cp[0].Accelerate(step,-maxspeed);

	if (keys[SDLK_PAGEDOWN])
	cp[2].Accelerate(step,maxspeed);
	if (keys[SDLK_PAGEUP])
	cp[2].Accelerate(-step,-maxspeed);
	}

	*/
	return 1;

}


void ThreeDPlot::TimerFunc(double t) {
	campos[0].Timer(t);
	campos[1].Timer(t);
	campos[2].Timer(t);
	campos[3].Timer(t);
	campos[4].Timer(t);
	campos[5].Timer(t);

}

void ThreeDPlot::Info() {

	cout << "[" << campos[0].value << " "
	     << campos[1].value << " "
	     << campos[2].value << "] "
	     << " phi=" << campos[3].value
	     << " psi=" << campos[4].value
	     << " fov"<< campos[5].value << endl;

}


void ThreeDPlot::DrawText00(double CurTime, double PlTime) const {


	//cout << "CCC " << CurTime << " " << PlTime << endl;

	if (CurTime>PlTime+1) {

		Point<int> pos = GetPos();

		glPushMatrix();
		glTranslatef(pos.x,pos.y,0.f);

		int width = GetWidth();
		int height = GetHeight();


		double a = (CurTime-PlTime-1);

		int rborder = 0;
		int lborder = 0;
		int tborder = 0;
		int bborder = 0;

		glBegin(GL_QUADS);

		//int z1 = (width-rborder-lborder)*.2;
		//int z2 = (height-tborder-bborder)*.2;

		//int z = min(z1,z2);
		//z = min(z1,10);



		glColor4f(0.f,0.f,0.f,a);
		glVertex3f(lborder,tborder,0.f);
		glVertex3f(width-rborder,tborder,0.f);
		glColor4f(0.f,0.f,0.f,a);
		glVertex3f(width-rborder,height-bborder,0.f);
		glVertex3f(lborder,height-bborder,0.f);
		glEnd();

		glPopMatrix();

	}


	float a = GetHue(); // [.3 ; .5]
	
	if (a>0) {
		Point<int> pos = GetPos();

		glPushMatrix();
		glTranslatef(pos.x,pos.y,0.f);

		int width = GetWidth();
		int height = GetHeight();

		int rborder = 8;
		int lborder = 8;
		int tborder = 8;
		int bborder = height-tborder-14;

		glBegin(GL_QUADS);

		//int z1 = (width-rborder-lborder)*.2;
		//int z2 = (height-tborder-bborder)*.2;

		//int z = min(z1,z2);
		//z = min(z1,10);

		glColor4f(.1f,.1f,.1f,a*.6f);
		glVertex3f(lborder,tborder,0);
		glVertex3f(width-rborder,tborder,0);
		glColor4f(.1f,.1f,.1f,a*.8f);
		glVertex3f(width-rborder,height-bborder,0);
		glVertex3f(lborder,height-bborder,0);
		//glColor4f(.1,.1,.1,a*.6);
		//glVertex3f(lborder,tborder,0);


		rborder = 14;
		lborder = 14;
		tborder = 14;
		bborder = height-tborder-10+8;

		glColor4f(.8f,.8f,.8f,a*.8f);
		glVertex3f(lborder,tborder,.1f);
		glVertex3f(width-rborder,tborder,.1f);
		glColor4f(.8f,.8f,.8f,a);
		glVertex3f(width-rborder,height-bborder,.1f);
		glVertex3f(lborder,height-bborder,.1f);
		
		//glColor4f(.8,.8,.8,a*.8);
		//glVertex3f(lborder,tborder,.1);


		double ff = CurTime/PlTime;
		if (ff>1) ff=1.;
		
		double p = lborder + ff*(width-rborder-lborder);
		float www = 4;
		float hhh = 3;
		//cout << CurTime << " " << PlTime << " " <<p << endl;

		glColor4f(.8f,0.f,0.f,a*.8f);
		glVertex3f(p-www,tborder-hhh,.2f);
		glVertex3f(p+www,tborder-hhh,.2f);
		glColor4f(.8f,0.f,0.f,a);
		glVertex3f(p+www,height-bborder+hhh,.2f);
		glVertex3f(p-www,height-bborder+hhh,.2f);
		//glColor4f(.8,.0,.0,a);
		//glVertex3f(p-www,tborder,.2);


		glEnd();


		/*
			glBegin(GL_LINE_STRIP);
			glColor4f(0,0,0,.5);
			glVertex3f(lborder,tborder,0.1);
			glVertex3f(width-rborder,tborder,0.1);
			glVertex3f(width-rborder,height-bborder,0.1);
			glVertex3f(lborder,height-bborder,0.1);
			glVertex3f(lborder,tborder,0.1);
			glEnd();

			glBegin(GL_LINE_STRIP);
			glColor4f(0,0,0,.1);
			glVertex3f(0,0,0.1);
			glVertex3f(0,height,0.1);
			glVertex3f(width,height,0.1);
			glVertex3f(width,0,0.1);
			glVertex3f(0,0,0.1);
			glEnd();*/

		glPopMatrix();
	}
}


void ThreeDPlot::drawFrame() const {

	int width = GetWidth();
	int height = GetHeight();
/*
	glBegin(GL_QUADS);
	float a = GetHue()*.2+.3; // [.3 ; .5]

	int z1 = (width-RBorder-LBorder)*.2;
	int z2 = (height-TBorder-BBorder)*.2;

	int z = min(z1,z2);
	//z = min(z1,10);
	
	glColor4f(1,1,1,a);
	glVertex3f(LBorder,TBorder,0);
	glVertex3f(width-RBorder,TBorder,0);
	glColor4f(1,1,1,a+.4);
	glVertex3f(width-RBorder,height-BBorder,0);
	glVertex3f(LBorder,height-BBorder,0);
	glColor4f(1,1,1,a);
	glVertex3f(LBorder,TBorder,0);
	glEnd();*/


	double step = 25*(High2-Low2)/(width-RBorder-LBorder);
	step = roundStep(step);


	//float col[3][4] = {{0.f,0.f,0.f,.4f},{0.f,0.f,0.f,.7f},{0.f,0.f,0.f,.1f}};


	//SetForeground(grey);


	/*
	glBegin(GL_LINE_STRIP);
	glColor4f(0,0,0,.5);
	glVertex3f(LBorder,TBorder,0.1);
	glVertex3f(width-RBorder,TBorder,0.1);
	glVertex3f(width-RBorder,height-BBorder,0.1);
	glVertex3f(LBorder,height-BBorder,0.1);
	glVertex3f(LBorder,TBorder,0.1);
	glEnd();*/


	glBegin(GL_LINE_STRIP);
	glColor4f(0.f,0.f,0.f,.1f);
	glVertex3f(0.f,0.f,0.1f);
	glVertex3f(0.f,height,0.1f);
	glVertex3f(width,height,0.1f);
	glVertex3f(width,0.f,0.1f);
	glVertex3f(0.f,0.f,0.1f);
	glEnd();

}

void ThreeDPlot::epsData(DataStorage&, DataStorage&, EpsWriter*) {}

int ThreeDPlot::rawHighLow(DataStorage&) {
	return 0;
}

void ThreeDPlot::drawText() const {}

}
