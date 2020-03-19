#pragma once

#include "tool.h"

namespace tw {

/** @ingroup toolbase
 *  @brief Rechteck mit Text
 */
class ToolRect {
public:

	/** Constructor. */
	ToolRect(int x_, int y_, int width_, int height_=14)
			: x(x_), y(y_), width(width_), height(height_), border(0), padding(2) {}

	bool contains(const Point<int> &p) const;
	bool contains(int x_,int y_) const;
	void Draw(double b = .1, float addz = 0.) const;

	int x,y,width,height;
	std::string text;
	int border;
	int padding;
};

}
