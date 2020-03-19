//
// C++ Interface: toolbox
//
// Description:
//
//
// Author: Matthias Knauer,,, <tulio@visurgis>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef TOOLWINDOW_H
#define TOOLWINDOW_H

#include "tool.h"


/** @defgroup toolbase Elements for Graphical User Interface
 *  @brief Button, Listbox, Edit, ...
 */

namespace tw {


/** @ingroup toolbase
 *  @brief ToolIcon.
 *
 */
class ToolIcon {
public:
	enum action{OPEN,WRITETXT,WRITEPNG,CLOSE};

	/** Constructor. */
	ToolIcon(action id_, int image_);

	void Draw(int width,int pos, bool moving, bool flag);
	int tx, ty;
	action id;

};


/** @ingroup toolbase
 *  @brief ToolWindow.
 *
 */
class ToolWindow {
public:
	/** Constructor. */
	ToolWindow(const std::string &t);

	/** Destructor. */
	virtual ~ToolWindow();

	void InitPos(int align_, int open_);
	void move(int ww, int hh);

	virtual void Init();

	virtual void Draw(bool isActive, int HT, const Point<int> &mouse, bool border=true);

	virtual std::string GetFilename();

	void AddIcon(ToolIcon::action id);
	void IconAction(ToolIcon::action id);
	virtual void WriteText();
	virtual void WriteImage();

	bool MouseMotion(bool isActive, SDL_MouseMotionEvent &m);
	bool MouseButton(SDL_MouseButtonEvent &m);

	virtual bool KeyboardFunc(SDL_Keysym &Keysym);
	// virtual bool SpecialKeys(Uint8 *keys);

	void RenderWindow();

	virtual void Resize();

protected:
	virtual bool mousebutton(int mx, int my, int button, int type);
	virtual bool mousemotion(int mx, int my);

	void innerProj(int HT);

public:
	std::string title;
	int x,y;
	int width,height;
	int scrollmax,scrollcur,scrolling;

	int open;
	
	int resizing,minsizex,minsizey;
	int resizeable;

protected:
	int theight;
	int align;

	int moving;
	Point<int> mouse;
	int setActiveModeTo;

public:
	int closeme;

	static void Rearrange() {
		posleft = 21;
		posright = 21;
	}

	int visible;


private:
	static int posleft;
	static int posright;

	static int scrollht;
	static int scrollwd;

	std::vector<ToolIcon*> icons;

};

}


#endif
