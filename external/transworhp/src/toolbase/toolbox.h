#pragma once

#include "toolwindow.h"

#include "toolitem.h"
#include "../base/dynloader.h"

#include "toolbutton.h"

#include <vector>

namespace tw {

class ToolItem;

/** @ingroup toolbase
 *  @brief Toolbox.
 *
 */
class ToolBox : public ToolWindow {
public:
	/** Constructor. */
	ToolBox(const std::string &t="");

	/** Destructor. */
	virtual ~ToolBox();

	REGISTER_CLASS(ToolBox)

	void Draw(bool isActive, int HT, const Point<int> &mouse, bool border=true);
	bool KeyboardFunc(SDL_Keysym &keysym);
	//bool SpecialKeys(Uint8 *keys);

	void Apply();
	void Discard();
	void Init();

	virtual void Update() {}

	std::vector<ToolButton*> buttons;
	std::vector<ToolItem*> items;


protected:
	bool mousebutton(int mx, int my, int button, int type);
	bool mousemotion(int mx, int my);
	void addStateButton(int &s, int id_, int image_, int modes_=0);
	void addButton(int id_, int image_);
	virtual void drawItem() {}

	std::vector<ToolItem *>::const_iterator activeItem;
	bool hasactiveitem;
};

}
