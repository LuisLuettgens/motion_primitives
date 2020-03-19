//
// Author: Matthias Knauer <knauer@math.uni-bremen.de>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
#ifdef WIN32
#include "windows.h"
#endif
#include "vektor.h"

#include <cmath>

template <class T>
Vektor<T> Vektor<T>::operator+(const Vektor<T> &other) const {
	Vektor ret;
	ret.x = x+other.x;
	ret.y = y+other.y;
	ret.z = z+other.z;
	return ret;
}

template <class T>
Vektor<T>& Vektor<T>::operator+=(const Vektor<T> &other) {
	x += other.x;
	y += other.y;
	z += other.z;
	return *this;
}

template <class T>
Vektor<T> Vektor<T>::operator-(const Vektor<T> &other) const {
	Vektor ret;
	ret.x = x-other.x;
	ret.y = y-other.y;
	ret.z = z-other.z;
	return ret;
}

template <class T>
Vektor<T>& Vektor<T>::operator-=(const Vektor<T> &other) {
	x -= other.x;
	y -= other.y;
	z -= other.z;
	return *this;
}

template <class T>
T Vektor<T>::Length() const {
	return std::sqrt(x * x + y * y + z * z);
}

template <class T>
Vektor<T> Vektor<T>::operator*(T other) const {
	Vektor ret;
	ret.x = x*other;
	ret.y = y*other;
	ret.z = z*other;
	return ret;
}

template <class T>
T Vektor<T>::operator*(const Vektor<T> &other) const {
	return x*other.x + y*other.y + z*other.z;
}

// inner product
template <class T>
Vektor<T> Vektor<T>::operator|(const Vektor<T> &other) const {
	return Vektor<T>(x*other.x,y*other.y,z*other.z);
}

template <class T>
Vektor<T>& Vektor<T>::operator*=(T other) {
	x *= other;
	y *= other;
	z *= other;
	return *this;
}

template <class T>
Vektor<T> Vektor<T>::operator/(T other) const {
	Vektor ret;
	ret.x = x/other;
	ret.y = y/other;
	ret.z = z/other;
	return ret;
}

template <class T>
Vektor<T>& Vektor<T>::operator/=(T other) {
	x /= other;
	y /= other;
	z /= other;
	return *this;
}

template <class T>
Vektor<T> Vektor<T>::rotX(T angle) const {
	Vektor ret;
	ret.x = x;
	ret.y = y*cos(angle)-z*sin(angle);
	ret.z = y*sin(angle)+z*cos(angle);
	return ret;
}

template <class T>
Vektor<T> Vektor<T>::rotY(T angle) const {
	Vektor ret;
	ret.x = x*cos(angle)+z*sin(angle);
	ret.y = y;
	ret.z = -x*sin(angle)+z*cos(angle);
	return ret;
}

template <class T>
Vektor<T> Vektor<T>::rotZ(T angle) const {
	Vektor ret;
	ret.x = x*cos(angle)-y*sin(angle);
	ret.y = x*sin(angle)+y*cos(angle);
	ret.z = z;
	return ret;
}

template <class T>
Vektor<T> Vektor<T>::cross(const Vektor<T> &other) const {
	Vektor ret;
	ret.x = y*other.z - z*other.y;
	ret.y = z*other.x - x*other.z;
	ret.z = x*other.y - y*other.x;
	return ret;
}

template <class T>
void Vektor<T>::ReduceToUnit() {

	T length = Length();
	if (length == T(0)) { length = T(1); }
	*this /= length;

}


template <>
/**
 * Get vector orthogonal to given triangle.
 * @param u Corner of triangle.
 * @param v Corner of triangle.
 * @param w Corner of triangle.
 * @param out Vector orthogonal to u-v and v-w.
 */
void calcNormal(const Vektor<float> &u, const Vektor<float> &v, const Vektor<float> &w, Vektor<float> &out) {

	out = (u-v).cross(v-w);
	out.ReduceToUnit();

}


template <class T>
Vektor<T> Vektor<T>::rotate(const Vektor<T> &axis, T angle) const {

	T c = 1 / std::sqrt(axis.y * axis.y + axis.z * axis.z);

	Vektor a1(   1    *x,
	             axis.z*c*y - axis.y*c*z,
	             axis.y*c*y + axis.z*c*z);

	T d = 1/axis.Length();

	Vektor a2(1/(c*d) *a1.x - axis.x/d*a1.z,
	          a1.y,
	          1/(c*d) *a1.z + axis.x/d*a1.x );

	Vektor a3 = a2.rotZ(angle);

	Vektor a4(1/(c*d) *a3.x + axis.x/d*a3.z,
	          a3.y ,
	          1/(c*d) *a3.z - axis.x/d*a3.x );

	Vektor a5(   1    *a4.x,
	              axis.z*c*a4.y + axis.y*c*a4.z,
	              -axis.y*c*a4.y + axis.z*c*a4.z);

	return a5;
}


// template <class T>
// Vektor<T> abs(const Vektor<T> &u) {
// 	return Vektor<T>(abs(u.x), abs(u.y), abs(u.z));
// }
//
// template <class T>
// Vektor<T> log(const Vektor<T> &u) {
// 	return Vektor<T>(log(u.x), log(u.y), log(u.z));
// }


template class Vektor<double>;
template class Vektor<float>;
