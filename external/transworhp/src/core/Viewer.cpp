#include "Viewer.h"

#include "../base/defines.h"
#include "../base/exception.h"
#include "../base/vectortools.h"

#include "xmlio.h"

#include "../gui/sdlcursor.h"

#include "../glbase/glfont.h"
#include "../glbase/viewport.h"
#include "../glbase/globject.h"

#include "../glplot/generalplot.h"
#include "../glplot/sparseplot.h"
#include "../glplot/adjPlot.h"

#include "../toolbase/tool.h"
#include "../toolbase/toolbox.h"
#include "../toolbase/toolmenuentry.h"

#include "../tool/toolregister.h"

#include "TransWORHP.h"
#include "TWfolder.h"
#include "TWparameter.h"
#include "worhp_info.h"
#include "TWdebug.h"
#include "TWproblem.h"

#include "tex_unihb.h"
#include "tex_zetem.h"
#include "tex_steinbeis.h"

#include "worhp/worhp.h"

#include <SDL2/SDL.h>
#include <iomanip>



using namespace std;

double clamp(double val, double low, double upp) {

	if (val < low) return low;
	if (val > upp) return upp;
	return val;
}

namespace tw {

Viewer::Viewer(TWparameter *twparam_) :
	SDLFrame(new SDLScreen(&twparam_->twwindow)),
	mouseOver(nullptr),
	running(false),
	twparam(twparam_),twfolder(nullptr),
	fps(0.0), framecount(0),
	animate(0),
	sceneTime(1000), playerTime(0.0),
	xlasttime(0.0),
	modifyPlot(nullptr),
	headermode(0)
{

	currentTime = SDL_GetTicks();
	startTime = currentTime;

	font = new GLFont();

	XMLNode *n = twparam_->xml->GetFirstChild("PATH");
	if (n) {
		path = n->GetText();
	}

	n = twparam_->xml->GetFirstChild("MODEL");
	while (n) {
		auto glo = new glObject();

		glo->Init(n);
		glo->InitObject(path);
		objects.push_back(glo);

		n = twparam_->xml->GetNextChild("MODEL");
	}

	registerTools();
	Tool::font = font;

	n = twparam_->xml->GetFirstChild("WINDOW");
	if (n) {
		XMLNode *n2 = n->GetFirstChild("MENU");
		if (n2) {
			Tool::InitStatic(n2);
		}

		XMLNode *nn = n->GetFirstChild("BACKGROUND");
		if (nn) {
			background.init(nn);
		}
		nn = n->GetFirstChild("LOGO");
		if (nn) {
			string s = nn->GetText();
			if (s == "ZeTeM") {
				background.logoID = 1;
			} else if (s == "UniBremen") {
				background.logoID = 2;
			} else if (s == "Steinbeis") {
				background.logoID = 3;
			} else {
				background.logoID = 0;
			}
		}
	}

	initTools();

	//hidemenu=1;
	menu.status.Clear();
	createMenu();

	SDLCursor::CreateCursors();
	background.initTexture();
}


Viewer::~Viewer() {

//	cout << "~Viewer" << endl;
	delete screen;
	delete font;

	plots.clear();

	SDLCursor::CleanCursors();

	for (DataStorage* ds : dsvector) {
		delete ds;
	}
	dsvector.clear();


#ifdef WIN32
	//cout << "~sdlquit" << endl;
#else
	SDL_Quit();
#endif
}





void Viewer::RenderScene() {

	screen->thethread->Lock();

	framecount++;

	int tex_width = screen->twwindow.size.x;
	int tex_height = screen->twwindow.size.y;

	///float mat_no_emission[] = {0.f, 0.f, 0.f, 0.0f};

	glColorMaterial(GL_FRONT, GL_DIFFUSE) ;
	glEnable(GL_COLOR_MATERIAL);
	//	 glMaterialf (GL_FRONT, GL_SHININESS, 100.0);
	//	 glMaterialf (GL_FRONT, GL_SPECULAR, 0,0,0,1);

	/*    glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient.GetData());
	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse.GetData());
	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular.GetData());
	glMaterialf(GL_FRONT, GL_SHININESS, shininess);
	glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, mat_emission.GetData());
	*/
	glEnable(GL_DEPTH_TEST);
	glDisable(GL_LINE_SMOOTH);

	glClearColor(background.color1.red(), background.color1.green(), background.color1.blue(), background.color1.alpha());
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


	glViewport(0, 0, tex_width, tex_height);

	//double nnear = 0.1;
	//double ffar = 30;

	glShadeModel(GL_SMOOTH);
	///(GL_BACK_LEFT);

	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);

	glEnable(GL_ALPHA_TEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glViewport(0, 0, screen->Width(), screen->Height());
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-screen->Width()/2, screen->Width()/2,
		-screen->Height()/2,screen->Height()/2,  -32.0, 10.0);
	glMatrixMode(GL_MODELVIEW);

	displayObjects();

	glClear(GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

#ifdef _MSC_VER
	//if (screen->twwindow.multisamplebuffers) glEnable(GL_MULTISAMPLE);
#endif
	display3dObjects();
#ifdef _MSC_VER
	//glDisable(GL_MULTISAMPLE);
#endif

	glClear(GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	displayText00();

	glDisable(GL_ALPHA_TEST);
	glDisable(GL_BLEND);

	glClear(GL_DEPTH_BUFFER_BIT);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	displayText();

	{
		glClear(GL_DEPTH_BUFFER_BIT);

		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();

		glOrtho(-screen->Width()/2, screen->Width()/2,
			-screen->Height()/2,screen->Height()/2,  -32.0, 10.0);
		glMatrixMode( GL_MODELVIEW );
		glLoadIdentity();

		glPushMatrix();

		for (auto it = plots.cbegin();it!=plots.end();++it) {
			(*it)->DrawText00(playerTime, sceneTime*.001);
		}
		glPopMatrix();
	}

	glDisable(GL_BLEND);
	glFlush();

	screen->thethread->Unlock();
}


void Viewer::displayObjects() {

	glPushMatrix();

	background.display(screen->Width(), screen->Height());
	background.displayLogo(screen->Width(), screen->Height());

	for (const auto &plot : plots) {

		plot->drawText();

		auto p = plot->GetPos();

		const GLint x = p.x + screen->Width()/2;
		const GLint y = p.y + screen->Height()/2;
		const GLint w = plot->GetWidth();
		const GLint h = plot->GetHeight();

		glScissor(x, y, w, h);
		glEnable(GL_SCISSOR_TEST);

		int ii = 0;
		for (DataStorage* data : dsvector) {
			plot->Draw(*data, dstop, ii);
			++ii;
			
		}
		glTranslatef(0.0f, 0.0f, 0.1f);

		glDisable(GL_SCISSOR_TEST);
	}

	background.displayMenu(screen->Width(),screen->Height());

	//optionsdlg->Draw(width,height);

	glPopMatrix();
}


void Viewer::display3dObjects() {

	if (x_now.empty()) {
		x_now.resize((*twfolder->phases.begin())->n_ode);
	}

	if (animate) {

		const double tt = playerTime/(sceneTime/1000.);

		(*twfolder->phases.begin())->solver->GetState(x_now.data(), tt);

		change=100;
		// if player aktive...

		//cout << playerTime << " " << sceneTime << endl;
	} else {
		const double tt = playerTime/(sceneTime/1000.);

		(*twfolder->phases.begin())->solver->GetState(x_now.data(), tt);
	}

	for (auto it = plots.cbegin(); it != plots.end(); ++it) {
		(*it)->DrawMore(this, x_now.data(), playerTime);
	}
}


void Viewer::displayText00() {

	glViewport(0, 0, screen->Width(), screen->Height());
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0, screen->Width(), -screen->Height(), 0, -1.0, 1.0);
	glMatrixMode(GL_MODELVIEW);

	if (headermode==0) {
		glColor4f(0.0f,0.0f,0.0f,0.8f);
	} else {
		glColor4f(1.0f,0.0f,0.0f,0.8f);
		background.color2 = color4(1.0f,0.0f,0.0f,0.3f);
	}

	const int ww = Tool::font->StringLength(header.c_str());
	font->printString(header.c_str(),static_cast<float>((screen->Width()-ww)/2),-16.0f,-0.5f);
}


void Viewer::displayText() {
	//glActiveTexture(GL_TEXTURE0);

	glClear(GL_DEPTH_BUFFER_BIT);
	glDisable(GL_LIGHTING);

	glViewport(screen->twwindow.reference.x, screen->twwindow.reference.y, screen->twwindow.size.x, screen->twwindow.size.y);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0.0, screen->twwindow.size.x, -screen->twwindow.size.y, 0.0, -5.0, 5.0);
	glMatrixMode(GL_MODELVIEW);

	glLineWidth(1.0f);

	glColor4f(1.0f,1.0f,1.0f,1.0f);

	glViewport(0, 0, screen->twwindow.size.x, screen->twwindow.size.y);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0.0, screen->twwindow.size.x, -screen->twwindow.size.y, 0.0, -12.0, 2.0);
	glMatrixMode(GL_MODELVIEW);

	glLineWidth(1.0f);
	glDisable(GL_LINE_SMOOTH);

	// FPS ausgeben
	std::string buf(" ");
	
	if (change > 1) {
		std::stringstream ss;
		ss << std::fixed << std::setprecision(2) << fps;
		buf = "fps: " + ss.str();
	}

	menu.infotext = std::move(buf);
	menu.Draw(screen->twwindow.size.x);
}


void Viewer::TimerFunc() {

	const int time1 = SDL_GetTicks();
	static int fpstime = time1;

	// FPS berechnen
	if (time1-fpstime > 500) {
		fps = framecount*1000.0/(time1-fpstime);
		fpstime = time1;
		framecount = 0;
	}

	if (xlasttime != 0.0) {

		const double t = (time1-xlasttime)/1000.0;

		if (animate) {
			playerTime += t;

			if (playerTime > 2 + sceneTime/1000.) {
				playerTime = 0.0;
			}
		}

		for (auto it = plots.cbegin(); it != plots.end(); ++it) {
			(*it)->TimerFunc(t);
		}
	}

	xlasttime = time1;

	currentTime = time1;
	for (size_t i = 0; i < numPlots(); i++) {
		getPlot(i)->Timer(this);
	}

	screen->thethread->Lock();
	menu.Timer(time1);
	screen->thethread->Unlock();

	if (currentTime-startTime > waittime * 1000 && !running) {
		running = true;
		screen->thethread->Run(twfolder, nullptr);
	}
}


void Viewer::writeZen() {

	static const std::vector<char> var = {'X','G','F','L','M'};
	static const std::vector<char> pert = {'P','R','Q','B'};

	for (char v : var) {
		//os << "DX/DP" << endl;
		ofstream os("ZenD" + std::string(1,v) + ".dat");

		for (char p : pert) {
			//os << "DX/DP" << endl;
			twfolder->WriteSensitivity(os, v, p);
		}
	}
}


void Viewer::on_SDL_QUIT() {
	twfolder->Interrupt = true;
}


void Viewer::Message(int id, const void */*param*/) {

	switch (id) {

		case 101:
			running = true;
			screen->thethread->Run(twfolder,nullptr);
			break;
		case 102:
			twfolder->Interrupt = true;
			break;

		case 103:
			//ph->ToMATLAB("transworhp.mat");
			break;

		case 104:
			{
				for (auto it = plots.cbegin(); it != plots.end(); ++it) {
					DataStorage *ds = *dsvector.begin();

					(*it)->Print(*ds,dstop);
				}
			}

			break;

		case 110:
			SDL_ShowCursor(SDL_DISABLE);
			MyQuit(-100);
			break;

		case 201:
			tilePlots();
			break;

		case 202:
			BasePlot::allplotnames ^= 1;
			change = 100;
			break;

		case 203:
			animate ^= 1;
			change = 100;
			break;

		case 204:
			{
				for (auto it = plots.cbegin(); it != plots.end(); ++it) {
					(*it)->Info();
				}
			}
			break;

		case 205:
			screen->twwindow.windowmode^=1;
			if (screen->twwindow.windowmode==0) change = 100;
			break;

		case 206:
			screen->twwindow.windowmode=10;
			change = 100;
			break;

		case 301:// NLP-Beschraenkungen
			TWdebug::PrintNLPconstraints(twfolder);
			break;

		case 302: // Opt-St-Zustaende
			(*twfolder->phases.begin())->solver->PrintOCPstates();
			break;

		case 401:
			TWdebug::PrintDF(twfolder,nullptr);
			break;

		case 402:
			TWdebug::PrintDG(twfolder,nullptr);
			break;

		case 403:
			TWdebug::PrintHM(twfolder,nullptr);
			break;

		case 404:  // MU, Lambda
			(*twfolder->phases.begin())->solver->PrintMultipliers();
			break;

		case 501:  // Zen
			writeZen();
			break;

		case 502: // Live Zen
		//	TWdebug::PrintHM(twfolder,0);
			break;

	}
}


bool Viewer::KeyboardFunc(SDL_Keysym &keysym) {

	if (menu.KeyboardFunc(keysym)) {
		return false;
	}

	return false;
}


bool Viewer::SpecialKeys(const Uint8 *keys) {

	if (menu.SpecialKeys(keys)) {
		return false;
	}

	/*
	SmoothMovement *cp = use_Cam? campos_close: campos;
	*/
	//int step = 3;
	//int maxspeed = 25;


	if (mouseOver) {
		if (mouseOver->SpecialKeys(keys)) {
			change = 100;
		}
	}

	return true;
}


bool Viewer::MouseButton(SDL_MouseButtonEvent &m) {

	if (menu.MouseButton(m)) {
		return false;
	}

	mouseFunc(m.button, m.type, Point<int>(m.x,m.y));

	return false;
}


void Viewer::MouseWheel(const SDL_MouseWheelEvent &m) {

	int x, y;
	SDL_GetMouseState(&x,&y);
	const Point<int> prel = Point<int>(x-screen->Width()/2,screen->Height()/2-y);

	if (mouseOver) {
		// realtive Koordinaten im Plot
		auto p = mouseOver->GetPos();
		const int x = p.x + mouseOver->LBorder;
		const int y = p.y;

		mouseOver->MouseWheelInput(m.y, Point<int>(prel.x-x, prel.y-y));
	}

	Redraw();
}


bool Viewer::MouseMotion(SDL_MouseMotionEvent &m) {

	if (menu.MouseMotion(m)) {
		return false;
	}

	if ((SDL_GetMouseState(nullptr, nullptr) & SDL_BUTTON_LMASK) ||
	    (SDL_GetMouseState(nullptr, nullptr) & SDL_BUTTON_RMASK)) {
		motionFunc(Point<int>(m.x,m.y));
	} else {
		passiveMotionFunc(Point<int>(m.x,m.y));
	}

	return false;
}


void Viewer::initTools() {

	Tool::font = font;
	Tool::frame = this;
}


void Viewer::createMenu() {

	// 100 = File
	ToolMenuEntry *m1 = menu.AddMenu("&File");

	menu.AddMenu(m1,"&Interrupt",102,"F3");
	menu.AddMenu(m1,"&Results to MATLAB File",103,"F4");

	menu.AddMenu(m1,"&Plot",104,"F1");

	/*menu.AddMenu(m1,"&Optimize",101,"F2");
	menu.AddSeparator(m1);*/
	menu.AddMenu(m1,"&Quit",110,"q");

	ToolMenuEntry *m2 = menu.AddMenu("&View");
	menu.AddMenu(m2,"&Tile Plots",201,"F5");
	menu.AddMenu(m2,"&Show Titles",202, "F6",&BasePlot::allplotnames);

	menu.AddSeparator(m2);
	menu.AddMenu(m2,"&Animation",203,"F7",&animate);
	menu.AddMenu(m2,"&Camera Info",204,"F8");

	menu.AddMenu(m2,"Freeze graphic output",205,"-", &screen->twwindow.windowmode);
	menu.AddMenu(m2,"Update graphic output",206,"#");


	ToolMenuEntry *m3 = menu.AddMenu("&Data");
	menu.AddMenu(m3,"&NLP Constraints",301);
	menu.AddMenu(m3,"&OCP States",302);
	/*	menu.AddMenu(m3,"D&F",303);
	menu.AddMenu(m3,"D&G",304);
	menu.AddMenu(m3,"&HM",305);*/
	menu.AddMenu(m3,"Multipliers",404,"F12");

	ToolMenuEntry *m4 = menu.AddMenu("Debu&g");
	menu.AddMenu(m4,"Check DF",401,"F9");
	menu.AddMenu(m4,"Check DG",402,"F10");
	menu.AddMenu(m4,"Check HM",403,"F11");

	ToolMenuEntry *m5 = menu.AddMenu("&Zen");
	menu.AddMenu(m5,"Calculate Zen Derivatives",501);
	menu.AddMenu(m5,"Live Zen", 502);
}


void Message(int id) {
	if (Tool::frame) {
		Tool::frame->Message(id);
	} else {
		//	gl->Message(id);
	}
}

void Message(int id, const void *param) {
	if (Tool::frame) {
		Tool::frame->Message(id,param);
	} else {
		//	gl->Message(id, param);
	}
}






// ---------------------------------------------------- //
//                  GLDisplay::Background
// ---------------------------------------------------- //

Viewer::Background::Background()
	:	logoID(1), flagGrid(false), waiting(false),
		mouseOverButton(0),
		border(3), border2(32),
		color1(0.7686f, 0.8824f, 1.0f, 1.0f),
		color2(0.6f,0.6f,0.7f,1.0f),
		colorGrid(0.4f,0.5f,0.6f,0.2f) {

}


void Viewer::Background::initTexture() {

	logotex[0].LoadTexture(zetem,0,128,128,4);
	logotex[1].LoadTexture(unihb,0,256,128,4);
	logotex[2].LoadTexture(steinbeis,0,256,128,4);
	//	tex[1].LoadTexture(menu,0,128,128,4);
}


void Viewer::Background::init(XMLNode *xml) {

	if (xml) {
		string s = xml->GetAttribute("color_bottom");
		if (!s.empty()) {
			color1 = s;
		}

		s = xml->GetAttribute("color_top");
		if (!s.empty()) {
			color2 = s;
		}

		s = xml->GetAttribute("color_grid");
		if (!s.empty()) {
			colorGrid = s;
			flagGrid = true;
		} else {
			flagGrid = false;
		}
	}
}


void Viewer::Background::displayMenu(int /*width*/, int /*height*/) {

}


void Viewer::Background::display(int width, int height) {

	glBegin(GL_QUADS);
	glColor4fv(color1.GetData());
	glVertex3f(-width/2.f+border,-height/2.f+border,-2.f);
	glVertex3f(width/2.f-border,-height/2.f+border,-2.f);
	glColor4fv(color2.GetData());
	glVertex3f(width/2.f-border,height/2.f-border2,-2.f);
	glVertex3f(-width/2.f+border,height/2.f-border2,-2.f);
	glColor4fv(color1.GetData());
	glVertex3f(-width/2.f+border,-height/2.f+border,-2.f);
	glEnd();

	if (flagGrid) {
		glBegin(GL_LINES);

		glColor4fv(colorGrid.GetData());
		for (int i = 0; i < width; i+=100) {
			glVertex3f(i,-height/2.f,-1.5f);
			glVertex3f(i,height/2.f,-1.5f);

			glVertex3f(-i,-height/2.f,-1.5f);
			glVertex3f(-i,height/2.f,-1.5f);
		}
		for (int i = 0; i < height; i+=100) {
			glVertex3f(-width/2.f,i,-1.5f);
			glVertex3f(width/2.f,i,-1.5f);
			glVertex3f(-width/2.f,-i,-1.5f);
			glVertex3f(width/2.f,-i,-1.5f);
		}
		glEnd();
	}
}


void Viewer::Background::displayLogo(int width, int height) {
	// Logo unten rechts

	if (logoID > 0) {
		glColor4f(1.f,1.f,1.f,.25f);
		logotex[logoID-1].Bind();
		glEnable( GL_TEXTURE_2D );
		glBegin(GL_QUADS);

		glNormal3f(0.f,0.f,1.f);

		float fx1 = 0.f;
		float fx2 = 1.f;
		float fy1 = 1.f;
		float fy2 = 0.f;

		int xx = logotex[logoID-1].GetWidth();
		int yy = logotex[logoID-1].GetHeight();

		glTexCoord2d(fx2,fy1);
		glVertex3f(width/2.f,-(height/2.f),-1.f);
		glTexCoord2d(fx2,fy2);
		glVertex3f(width/2.f,-(height/2.f-yy),-1.f);
		glTexCoord2d(fx1,fy2);
		glVertex3f(width/2.f-xx,-(height/2.f-yy),-1.f);
		glTexCoord2d(fx1,fy1);
		glVertex3f(width/2.f-xx,-(height/2.f),-1.f);

		glEnd();
		glDisable(GL_TEXTURE_2D);
	}
}

/////////////////// DragInfo ///////////////////

Viewer::DragInfo::DragInfo() : plot(nullptr) {}


void Viewer::DragInfo::set(BasePlot *p,int m, const Point<int> &pt) {
	plot = p;
	point = pt;
	mode = m;
	geom = p->geom.cur;
}


void Viewer::closeAll() {

	if (!twfolder) {
		return;
	}

	plots.clear();
}

void Viewer::update() {

	if (dsvector.empty()) {
		dsvector.reserve(twfolder->phases.size());
		for (size_t i = 0; i < twfolder->phases.size(); ++i) {
			auto ds = new DataStorage;
			dsvector.push_back(ds);
		}
	}

	auto it = twfolder->phases.begin();
	auto it2 = dsvector.begin();

	double t0 = +1e200;
	double tf = -1e200;

	for (int i = 0; it != twfolder->phases.end(); ++it, ++it2, i++) {

		update(*it,i,*it2);

		const double tf0 = (*it2)->getTF();
		if (tf0 > tf) {
			tf = tf0;
		}

		if ((*it)->solver->T[0] < t0) {
			t0 = (*it)->solver->T[0];
		}
	}

	change = 100;

	//ds.SetData(len,ndgl,nsteuer,nneben,t,x,u,g);

	for (const auto &plot : plots) {
		plot->SetMinTime(DataStorage::timemode,t0);
		plot->SetMaxTime(DataStorage::timemode,tf);
	}

	DataStorage::MaxTF = tf;
	DataStorage::MinT0 = t0;

	autoScale();
}


void Viewer::update(TransWorhpProblem *ph, int i, DataStorage *ds) {
	//	if (!ph) return;

	if (ph->freetime) {
		const double t = ph->solver->p(0);
		const double t0 = ph->t0;
		setFloatTime(t, t0);
	}

	ph->solver->updateViewer(ds, temptime[i]);
}


void Viewer::Matrix(string s, WorhpMatrix *m) {

	SparsePlot *plot = new SparsePlot(m);
	plot->SetEpsTitle(std::move(s));
	addPlot(plot);

	tilePlots();
}

void Viewer::init(TWfolder *p) {

	twfolder = p;

	TransWorhpProblem* ph = *twfolder->phases.begin();

	ph->selectWindows(this);

	setTimeIsTime();
	//noautochange_();

	temptime.resize(twfolder->phases.size());
	setTempTime();

	//  plot_(ph.n_dis, ph.n_ode, ph.n_ctrl, NNEBEN,ph.T,ph.X,ph.U,ph.X);


	string name("TransWORHP - " + ph->id);

	for (TransWorhpProblem* ph : twfolder->phases) {
		name += ph->solver->setViewerTitle();
	}

	SDL_SetWindowTitle(screen->window, name.c_str());

	//wait_();

	ph->OpenWindows(this);

	tilePlots();
}

void Viewer::AddStateView(const int n, string name) {

	GeneralPlot* plot = new GeneralPlot('x',n,static_cast<int>(numPlots()));

	if (name.empty()) {
		name = "X" + std::to_string(n);
	}
	plot->SetEpsTitle(std::move(name));

	for (TransWorhpProblem* ph : twfolder->phases) {
		// fixe Anfangs- und Endwerte
		vector<int> indices;
		indices.reserve(ph->n_ode*2); //Speicherplatz fuer AW und EW
		ph->solver->GetBoundaryIndices(indices, n);
		plot->dot_indices.push_back(move(indices));
	}

	addPlot(plot);
}

void Viewer::AddControlView(const int n, string name) {

	GeneralPlot *plot = new GeneralPlot('u',n,static_cast<int>(numPlots()));

	if (name.empty()) {
		name = string("U") + std::to_string(n);
	}
	plot->SetEpsTitle(std::move(name));

	plot->redflag = 1;
	//Redflag();

	addPlot(plot);
}

void Viewer::AddIntegralView(const int n, string name) {

	GeneralPlot *plot = new GeneralPlot('l', n, static_cast<int>(numPlots()));

	if (name.empty()) {
		name = string("Integral_") + std::to_string(n);
	}

	plot->SetEpsTitle(std::move(name));
	addPlot(plot);
}

void Viewer::AddLambdaView(const int n, string name) {

	TransWorhpProblem* ph = *twfolder->phases.begin();

	int m, type;

	if (n < ph->n_ode) {// Zustand
		m = n;
		type = 0;
	} else { // Steuerung
		m = n - ph->n_ode;
		type = 1;
	}

	string buf("Beschraenkungen " + std::move(name));

	lambdaPlot2* plot = new lambdaPlot2(buf, m, type, ph->solver->T.data(), ph->solver->Lambda, ph->n_dis, ph->n_ctrl, ph->n_ode, ph->solver->twdiscretization);
	plot->SetEpsTitle(std::move(buf));
	addPlot(plot);
}

void Viewer::AddMuView(const int n, string name) {

	TransWorhpProblem* ph = *twfolder->phases.begin();

	string buf("Adjungierte x" + std::to_string(n) + " " + std::move(name));

	adjPlot* plot = new adjPlot(n, ph->solver->T.data(), ph->solver->Mu, ph->n_dis, ph->n_ode, ph->solver->twdiscretization);
	plot->SetEpsTitle(std::move(buf));
	addPlot(plot);
}


void Viewer::selectWindows() {

	TransWorhpProblem* ph = *twfolder->phases.begin();

	for (int d = 0; d < ph->n_ode; d++) {
		string s( ph->GetXTitle(d) );
		if (s.empty()) {
			s = string("X") + std::to_string(d);
		}
		AddStateView(d, std::move(s));
	}

	for (int d = 0; d < ph->n_ctrl; d++) {
		string s( ph->GetUTitle(d) );
		if (s.empty()) {
			s = string("U") + std::to_string(d);
		}
		AddControlView(d, std::move(s));
	}

	for (int d = 0; d < ph->n_integral; d++) {
		string s ( "Integral" + std::to_string(d) );
		AddIntegralView(d, std::move(s));
	}
}


void Viewer::setTempTime() {

	int j = 0;

	for (TransWorhpProblem* ph : twfolder->phases) {

		ph->solver->setTemptimeForViewer(temptime[j]);
		j++;
	}
	update();
}


void Viewer::tilePlots(int TT) {

	size_t size = plots.size();
	int cols,rows;

	if (size <= 4) {
		cols = 2;
		rows = 2;
	} else if (size <= 6) {
		cols = 3;
		rows = 2;
	} else if (size <= 9) {
		cols = 3;
		rows = 3;
	} else if (size <= 12) {
		cols = 4;
		rows = 3;
	} else if (size <= 16) {
		cols = 4;
		rows = 4;
	} else if (size <= 20) {
		cols = 5;
		rows = 4;
	} else if (size <= 25) {
		cols = 5;
		rows = 5;
	} else if (size <= 30) {
		cols = 6;
		rows = 5;
	} else if (size <= 36) {
		cols = 6;
		rows = 6;
	} else { /*if (size<=56)*/
		cols = 8;
		rows = (size-1)/8 + 1;
	}

	const int newwidth = (screen->Width()-2*background.border)/cols;
	const int newheight = (screen->Height()-background.border-background.border2)/rows;
	//cout << width << " " << height << " " << newwidth << " " << newheight << endl;

	for (size_t i = 0; i < size; i++) {
		int x = static_cast<int>(i%cols)*(newwidth);
		int y = static_cast<int>(i/cols)*(newheight);

		//	cout << i << " " << (i%cols) << " " << (i/cols) << " " << x << " " << y << endl;
		//		plots[i]->newHighLow();
		//		plots[i]->calcYscale();

		BasePlot::Geometry a(newwidth*plots[i]->extx,newheight*plots[i]->exty);
		a.pos = Point<int>(x-screen->Width()/2+background.border,screen->Height()/2-y-newheight-background.border2);

		plots[i]->MoveTo(currentTime, a, TT);
	}

	autoScale();
}


void Viewer::autoScale() {

	for (size_t i = 0; i < numPlots(); i++) {

		double Lo_ = plots[i]->Low;
		double Lo2_ = plots[i]->Low2;
		double Hi_ = plots[i]->High;
		double Hi2_ = plots[i]->High2;

		for (auto it = dsvector.begin(); it != dsvector.end(); ++it) {

			plots[i]->newHighLow(**it);

			if (it == dsvector.begin()) {
				Lo_ = plots[i]->Low;
				Lo2_ = plots[i]->Low2;
				Hi_ = plots[i]->High;
				Hi2_ = plots[i]->High2;
			}
			else {
				if (plots[i]->Low<Lo_) Lo_ = plots[i]->Low;
				if (plots[i]->Low2<Lo2_) Lo2_ = plots[i]->Low2;
				if (plots[i]->High>Hi_) Hi_ = plots[i]->High;
				if (plots[i]->High2>Hi2_) Hi2_ = plots[i]->High2;
			}
		}

		plots[i]->Low = Lo_;
		plots[i]->Low2 = Lo2_;
		plots[i]->High = Hi_;
		plots[i]->High2 = Hi2_;

		plots[i]->calcYscale();
	}
}


void Viewer::addPlot(BasePlot* p) {
	plots.push_back(p);
}


size_t Viewer::numPlots() const {
	return plots.size();
}


BasePlot* Viewer::lastPlot() const {
	return plots.back();
}


BasePlot* Viewer::getPlot(int i) const {
	return plots[i];
}


//extern TransWorhp *globalTransWorhp;
extern TWfolder *globaltwfolder;

void* calcTransWorhp(int *) {

	if (globaltwfolder) globaltwfolder->WorhpLoop();
//	else if (globalTransWorhp) globalTransWorhp->WorhpLoop();

	return nullptr;
}


void Viewer::motionFunc(const Point<int> &p) {

	if (dragInfo.plot) {

		int dwidth = dragInfo.plot->GetWidth();
		int dheight = dragInfo.plot->GetHeight();
		Point<int> dpos = dragInfo.plot->GetPos();

		if (dragInfo.mode==1) {
			dpos = Point<int>(p.x,screen->Height()/2-p.y)-dragInfo.point+dragInfo.geom.pos;

			dpos.x = clamp(dpos.x,background.border-screen->Width()/2, screen->Width()/2-dwidth-background.border);
			dpos.y = clamp(dpos.y,background.border-screen->Height()/2, screen->Height()/2-dheight-background.border2);
		}
		if (dragInfo.mode&2) { // L
			int w = -(p.x-dragInfo.point.x)+dragInfo.geom.width;

			w = clamp(w,BasePlot::minwidth,screen->Width()/2+dragInfo.geom.pos.x+dragInfo.geom.width);
			dpos.x = dragInfo.geom.pos.x+dragInfo.geom.width-w;
			dwidth = w;
		}
		if (dragInfo.mode&4) { // R
			int w = p.x-dragInfo.point.x+dragInfo.geom.width;
			w = clamp(w,BasePlot::minwidth,screen->Width()/2-dragInfo.geom.pos.x);
			dwidth = w;
		}
		if (dragInfo.mode&8) { // T
			int h = screen->Height()/2-p.y-dragInfo.point.y+dragInfo.geom.height;
			h = clamp(h,BasePlot::minheight,screen->Height()/2-dragInfo.geom.pos.y);
			dheight =  h;
		}
		if (dragInfo.mode&16) { // B
			int h = -(screen->Height()/2-p.y-dragInfo.point.y)+dragInfo.geom.height;
			h = clamp(h,BasePlot::minheight,screen->Height()/2+dragInfo.geom.pos.y+dragInfo.geom.height);
			dpos.y = dragInfo.geom.pos.y+dragInfo.geom.height-h;
			dheight = h;
		}
		if (dragInfo.mode&(16+8+4+2+1)) {
			BasePlot::Geometry g(dwidth, dheight);
			g.pos = dpos;
			dragInfo.plot->MoveTo(currentTime, g, 0);
		}
		if (dragInfo.mode&(16+8+4+2)) {
			dragInfo.plot->calcYscale();
		}

		change = 100;

	}


	if (modifyPlot) {

		Point<int> prel = Point<int>(p.x-screen->Width()/2,screen->Height()/2-p.y);

		if (modifyPlot->containsArea(prel)) {
			DataStorage *ds = *dsvector.begin();

			if (modifyPlot->MouseInput(*ds,0,prel)) {

			}
		}

		change = 100;
	}
}

void Viewer::passiveMotionFunc(const Point<int> &p) {

	change = 100;

	mouseOver = nullptr;
	background.mouseOverButton = 0;

	int a = 0;

	for (auto it = plots.rbegin(); it != plots.rend(); ++it) {
		if ((*it)->contains(Point<int>(p.x-screen->Width()/2,screen->Height()/2-p.y))) {
			mouseOver = (*it);
			if (a!=-1)
				a = (*it)->containsArea(Point<int>(p.x-screen->Width()/2,screen->Height()/2-p.y));
			break;
		}
	}

	if (a==-1) a = 1;
	cursor[a].Set();
}

/*
1: Mitte
2: L
4: R
8: T
16: B
*/



void Viewer::mouseFunc(int button, int state, const Point<int> &p) {

	static int doppeltime = SDL_GetTicks();
	// damit erster Klick nicht als Doppelklick erkannt wird
	static bool first = true;

	Point<int> prel = Point<int>(p.x-screen->Width()/2,screen->Height()/2-p.y);

	if (modifyPlot && state == SDL_MOUSEBUTTONUP && button ==  SDL_BUTTON_LEFT) {
		DataStorage *ds = *dsvector.begin();

		modifyPlot->MouseInput(*ds ,2,prel);
	}

	dragInfo.plot = nullptr;
	modifyPlot = nullptr;

	// Doppelklick registrieren
	if (!first && state == SDL_MOUSEBUTTONDOWN) {
		int dtime = SDL_GetTicks();
		if (dtime-doppeltime < 250) {
			//std::cout << "DOPPELKLCIK" <<dtime-doppeltime << endl;

			for (auto it = plots.rbegin(); it != plots.rend(); ++it) {
				const int a =(*it)->containsArea(prel);
				if (a) {
					maxPlot(*it);//	dragInfo.Set(*it,a,Point(p.x,height/2-p.y));
					break;
				}
			}
		}
		doppeltime = dtime;
	}
	first = false;

	// Fenster verschieben
	if (button == SDL_BUTTON_RIGHT && state == SDL_MOUSEBUTTONDOWN) {

		for (auto it = plots.rbegin(); it != plots.rend(); ++it) {
			const int a =(*it)->containsArea(prel);
			if (a) {
				dragInfo.set(*it, a, Point<int>(p.x,screen->Height()/2-p.y));
				break;
			}
		}
	}

	if (button == SDL_BUTTON_LEFT && state == SDL_MOUSEBUTTONDOWN) {

		for (auto it = plots.rbegin(); it != plots.rend(); ++it) {
			const int a = (*it)->containsArea(prel);
			if (a == 1) {
				DataStorage *ds = *dsvector.begin();

				if ((*it)->MouseInput(*ds, 1, prel)) {
					modifyPlot = *it;
					break;
				}

			} else if (a) {
				dragInfo.set(*it, a, Point<int>(p.x, screen->Height()/2-p.y));
				break;
			}
		}
	}

	if (button == SDL_BUTTON_MIDDLE && state == SDL_MOUSEBUTTONUP) {
		tilePlots();
	}

	change = 100;
}

void Viewer::maxPlot(BasePlot *p) {

	int size = static_cast<int>(plots.size());
	int rows = (screen->Height()-background.border-background.border2)/BasePlot::minheight;
	int cols = size/rows+1;

	//	cout << "MAXPLOT" << rows << " " << cols<< endl;
	int newwidth = (screen->Width()-2*background.border)-BasePlot::minwidth*cols;
	int newheight = (screen->Height()-background.border-background.border2);

	BasePlot::Geometry a(newwidth,newheight);
	int x = 0;
	int y = 0;
	a.pos = Point<int>(x-screen->Width()/2+background.border,screen->Height()/2-y-newheight-background.border2);

	p->MoveTo(currentTime, a, 500);

	//cout << width << " " << height << " " << newwidth << " " << newheight << endl;

	for (int i = 0; i < size; i++) {
		if (plots[i] != p) {
			int x = newwidth + (i/rows)*(BasePlot::minwidth);
			int y = (i%rows)*(BasePlot::minheight);

			//	cout << i << " " << (i%cols) << " " << (i/cols) << " " << x << " " << y << endl;
			//	plots[i]->newHighLow();
			//	plots[i]->calcYscale();

			BasePlot::Geometry a(BasePlot::minwidth,BasePlot::minheight);
			a.pos = Point<int>(x-screen->Width()/2+background.border,screen->Height()/2-y-BasePlot::minheight-background.border2);

			plots[i]->MoveTo(currentTime, a, 500);
		}
	}
}


int Viewer::Redraw() {

	static int lastrender = 0;

	int ret = 0;

	if (screen->twwindow.windowmode==1) return ret;

	if (change) {
		if (SDL_GetTicks() - lastrender > 20) {
			change--;
			lastrender = SDL_GetTicks();
			ret=1;
		}
	}

	if (screen->twwindow.windowmode > 1) {
		screen->twwindow.windowmode--;
	}

	//	cout << "Redraw? " << ret << " " << change << endl;
	return ret;

}

void Viewer::PhasePlot(std::string s, Funktionenzeiger2 func, int d1, int d2) {

	if (!twfolder) return;

	Data2Plot* plot = new Data2Plot(s, func, d1, d2, static_cast<int>(numPlots()));

	plot->SetEpsTitle(std::move(s));
	addPlot(plot);
}


void Viewer::ThreeD(std::string s, XMLNode *xml, plot3d f) {

	ThreeDPlot* plot = new ThreeDPlot(xml, f);
	plot->SetEpsTitle(std::move(s));
	addPlot(plot);
}


void Viewer::Data(std::string s, Funktionenzeiger2 func, int index) {

	if (!twfolder) return;

	DataPlot* plot = new DataPlot(s, func, index, static_cast<int>(numPlots()));

	plot->SetEpsTitle(std::move(s));
	addPlot(plot);
}


void Viewer::AddDynCompareCurve(int i, int *cmpstep, double* cmp, double* cmptime, int *n) {

	if (!twfolder) return;

	plots[i]->AddDynCompareCurve(cmptime,cmp,cmpstep,n);
}


void Viewer::AddCompareCurve(int &cmpstep, double* cmp, double* cmptime, int n) {

	if (!twfolder) return;

	BasePlot *plot = lastPlot();
	if (plot) {
		plot->AddCompareCurve2(cmptime, cmp, cmpstep, n);
	}
}


void Viewer::AddCompareCurve(int i, double* cmptime, double* cmp, int &cmpstep, int n) {
	if (!twfolder) return;
	plots[i]->AddCompareCurve2(cmptime, cmp, cmpstep, n);
}


void Viewer::setHeader(char *s, int mode) {

	header = s;
	headermode = mode;
}

void Viewer::setFloatTime(double t, double t0) {

	DataStorage::timemode = DataStorage::fixed_float_as_time;
	DataStorage::floattime=t;
	DataStorage::floatstart=t0;

	sceneTime = static_cast<int>(t*1000);
}

void Viewer::setTimeIsTime() {

	DataStorage::timemode = DataStorage::time_is_time;
}


void Viewer::disFehler(const vector<vector<double>> &zeit, const vector<vector<double>> &fehler, const TWdiscretization *twdiscretization) {

	if (!twfolder) return;

	if (twparam->meshref_PLOT_SW) {
		std::string buf("Schrittweite");
		gitterPlot *plot = new gitterPlot(buf, 0, zeit, fehler);
		plot->SetEpsTitle(std::move(buf));
		addPlot(plot);
	}

	if (twparam->meshref_PLOT_ERR) {
		std::string buf("Diskretisierungsfehler");
		gitterPlot *plot = new gitterPlot(buf, 1, zeit, fehler);
		plot->SetEpsTitle(std::move(buf));
		addPlot(plot);
	}

	if (twparam->meshref_PLOT_MP) {
		std::string buf("Stuetzstellen");
		punktPlot *plot = new punktPlot(buf, 2, zeit, twdiscretization);
		plot->SetEpsTitle(std::move(buf));
		addPlot(plot);
	}

	tilePlots();
}

void Viewer::restrictionPlot(const double *zeit, const double *lambda, int n_dis, int n_ctrl, int n_ode, const TWdiscretization *twdiscretization) {

	if (!twfolder) return;

	if (twparam->meshref_PLOT_LAMBDA) {
		std::string buf("Beschraenkungen");
		lambdaPlot *plot = new lambdaPlot(buf, zeit, lambda, n_dis, n_ctrl, n_ode, twdiscretization);
		plot->SetEpsTitle(std::move(buf));
		addPlot(plot);
	}

	tilePlots();
}

}
