#pragma once

#include <GL/glew.h>
#include <GL/gl.h>

#ifndef _WIN32
#include "TWGUIconfig.h"
#endif

#include "../base/defines.h"
#include "../base/color4.h"

#include "../glbase/texture.h"

#include "../gui/sdlframe.h"
#include "../gui/sdlthread.h"

#include "../toolbase/toolmenu.h"

#include "../glplot/baseplot.h"
#include "../glplot/plot.h"
#include "../glplot/threedplot.h"
#include "../glplot/xopt_data.h"

#include <vector>

namespace tw {
	class glObject;
	class ToolBox;
	class TWfolder;
	class TransWorhpProblem;
}
class FFont;

#ifndef DllExport
#ifdef _WIN32
#define DllExport __declspec(dllexport)
#else
#define DllExport
#endif
#endif

#ifdef _MSC_VER
#pragma warning(disable : 4251)
#endif

#ifdef TRANSWORHP_GRAPHICS

namespace tw {

class DllExport Viewer : public SDLFrame {

public:
	Viewer(TWparameter *twparam);
	~Viewer();

	void RenderScene() override;

	void TimerFunc() override;
	bool SpecialKeys(const Uint8 *keys) override;
	bool KeyboardFunc(SDL_Keysym &keysym) override;

	bool MouseButton(SDL_MouseButtonEvent &m) override;
	void MouseWheel(const SDL_MouseWheelEvent &m) override;
	bool MouseMotion(SDL_MouseMotionEvent &m) override;

	void Message(int id, const void *param=nullptr) override;

	int Redraw() override;

	void on_SDL_QUIT() override;

	void init(TWfolder *p);

	/** Update graphics */
	void update();
	/** Close windows */
	void closeAll();

	void tilePlots(int TT=500);
	void autoScale();

	void selectWindows();

	void AddStateView(int n, std::string name);
	void AddControlView(int n, std::string name);
	void AddIntegralView(int n, std::string name);
	void AddLambdaView(int n, std::string name);
	void AddMuView(int n, std::string name);

	void AddCompareCurve(int &cmpstep, double* cmp, double* cmptime, int n);
	void AddDynCompareCurve(int i, int *cmpstep, double* cmp, double* cmptime, int *n);
	void AddCompareCurve(int i, double* cmptime, double* cmp, int &cmpstep, int n);

	void PhasePlot(std::string s, Funktionenzeiger2 func, int d1, int d2);
	void Data(std::string s, Funktionenzeiger2 func, int index);
	void ThreeD(std::string s, XMLNode *xml, plot3d f);
	void Matrix(std::string s, WorhpMatrix *m);

	/** number of plots */
	size_t numPlots() const;
	BasePlot* getPlot(int i) const;

	void writeZen();

	void setHeader(char *s, int mode);

	void disFehler(const std::vector<std::vector<double>> &zeit, const std::vector<std::vector<double>> &fehler, const TWdiscretization *twdiscretization);
	void restrictionPlot(const double *zeit, const double *lambda, int n_dis, int n_ctrl, int n_ode, const TWdiscretization *twdiscretization);

private:

	void initTools();
	void createMenu();

	void mouseFunc(int button, int state, const Point<int> &p);
	void passiveMotionFunc(const Point<int> &p);
	void motionFunc(const Point<int> &p);

	void setFloatTime(double t, double t0);
	void setTimeIsTime();

	void update(TransWorhpProblem *ph, int i, DataStorage *ds);

	void setTempTime();

	void maxPlot(BasePlot *p);
	void addPlot(BasePlot *p);

	BasePlot* lastPlot() const;

	void displayObjects();
	void display3dObjects();

	void displayText();
	void displayText00();

public:
	BasePlot *mouseOver;

	int startTime;
	int currentTime;
	int change;
	bool running;

private:

	// std::unique_ptr macht hier Probleme unter Visual Studio
	std::vector<BasePlot *> plots;

	TWparameter *twparam;
	TWfolder *twfolder;

	FFont *font;

	double fps;
	int framecount;
	int animate;
	int sceneTime;
	double playerTime;

	std::string path;

	double xlasttime;

	std::vector<DataStorage*> dsvector;
	DataStorage dstop;

	std::vector<std::vector<double>> temptime;

	BasePlot *modifyPlot;

	std::vector<glObject*> objects;

	ToolMenu menu;

	std::vector<double> x_now;

	std::string header;
	int headermode;

	struct Background {
		Background();
		void init(XMLNode *xml);
		void initTexture();
		void displayMenu(int width, int height);
		void displayLogo(int width, int height);
		void display(int width, int height);

		// Logo unten rechts
		int logoID;
		// Gitternetz im Hintergrund
		bool flagGrid;

		bool waiting;
		int mouseOverButton;
		int border, border2;
		Texture logotex[3];
		color4 color1;
		color4 color2;
		color4 colorGrid;
	}
	background;

	struct DragInfo {

		DragInfo();
		void set(BasePlot *p,int m, const Point<int> &pt);

		BasePlot *plot;
		Point<int> point;
		int mode;
		BasePlot::Geometry geom;
	}
	dragInfo;
};

}

#endif
