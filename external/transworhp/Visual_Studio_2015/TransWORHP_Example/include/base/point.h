//
// Author: Matthias Knauer <knauer@math.uni-bremen.de>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
#pragma once

#include <iostream>
#include <cmath>

/** @defgroup glbase Basics for OpenGL
 *  @brief ...
 */

namespace tw {

/** @ingroup glbase
 *  @brief Point with (x/y) integer coordinates.
 * 
 */
template <typename T>
struct Point {

	Point();
	Point(T x_, T y_);

	Point operator+(const Point &other) const {
		Point p;
		p.x = x+other.x;
		p.y = y+other.y;
		return p;
	}
	Point operator*(T v) const {
		Point p;
		p.x = x*v;
		p.y = y*v;
		return p;
	}
	Point operator-(const Point &other) const {
		Point p;
		p.x = x-other.x;
		p.y = y-other.y;
		return p;
	}
	Point ReduceToUnit() const {
		T len = static_cast<T>(std::sqrt(static_cast<double>(x*x + y*y)));
		
		if (len>1e-6) return Point(x/len,y/len);
		return Point(0,0);
	}
	Point rotleft() const {
		return Point(-y,x);
	}
	Point rotright() const {
		return Point(y,-x);
	}
	
	bool isin(const Point &c1, const Point &c2) const;

	T *data() {
		return &x;
	}

	T x,y;


	friend std::ostream &operator<<(std::ostream &os, const Point<T> &p) {

		os << "[" << p.x << ";"<< p.y << "]";
		return os;
	}
	
	bool operator<(const Point<T> &other) const;

};

}
