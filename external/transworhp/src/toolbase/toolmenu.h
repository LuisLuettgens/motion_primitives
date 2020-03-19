#pragma once

#include "../base/point.h"
#include "SDL2/SDL.h"
#include "toolmenukey.h"
#include "toolstatus.h"
#include "../glbase/smoothmovement.h"

#include <vector>
#include <string>

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

namespace tw {

class ToolMenuEntry;

/** @ingroup glbase
 *  @brief Menu
 *
 */
class DllExport ToolMenu {

public:
	ToolMenu();
	virtual ~ToolMenu();
	void Draw(int w) const;

	ToolMenuEntry *AddMenu(const std::string &s, int id=0);
	//MenuEntry *AddMenu(MenuEntry *th, const std::string &s, int id);
	ToolMenuEntry *AddMenu(ToolMenuEntry *th, const std::string &s, int id, const std::string &hot="", int *selected=0);
	ToolMenuEntry *AddMenu(ToolMenuEntry *th, const std::string &s, int id, float *val, float minval, float maxval, const std::string &hot="");
	ToolMenuEntry *AddMenu(ToolMenuEntry *th, const std::string &s, int id, SmoothMovement &sm, const std::string &hot="");

	ToolMenuEntry *AddSeparator(ToolMenuEntry *th);

	bool MouseButton(SDL_MouseButtonEvent &m);
	bool MouseMotion(SDL_MouseMotionEvent &m);
	void Timer(int t);
	bool KeyboardFunc(SDL_Keysym &keysym);
	bool SpecialKeys(const Uint8 *keys);

	ToolStatus status;

	std::string infotext;

private:
	std::vector<ToolMenuEntry*> me;
	ToolMenuEntry *openentry;
	int opened;

	int visible;
	int yy;

	std::vector<ToolMenuKey> mk;

	int menuposct;

};

}
