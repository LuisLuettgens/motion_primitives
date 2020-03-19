//
// C++ Implementation: point
//
// Description:
//
//
// Author: Matthias Knauer <knauer@math.uni-bremen.de>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "point.h"

namespace tw {

template <class T>
Point<T>::Point() : x(0),y(0) {}

template <class T>
Point<T>::Point(T x_, T y_) : x(x_),y(y_) {}

template <class T>
bool Point<T>::isin(const Point<T> &c1, const Point<T> &c2) const {

	return x >= c1.x && x <= c2.x && y >= c1.y && y <= c2.y;
}

template <class T>
bool Point<T>::operator<(const Point<T> &other) const {
	if (x == other.x) { return y < other.y; }
	return x<other.x;
}

template struct Point<double>;
template struct Point<float>;
template struct Point<int>;

} // namespace tw
