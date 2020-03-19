//
// Author: Matthias Knauer <knauer@math.uni-bremen.de>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
#pragma once

#include <vector>
#include <iostream>
#include <set>

namespace tw {

#ifndef TYPESAFE

/** @name Access templates for STL classes.  */
// @{

/** @ingroup base
 *  Empty vector and clear pointers
 *  @param vlist
 */
template <class T>
void ClearPointer(std::vector<T *> &vlist) {
	typename std::vector<T *>::iterator it = vlist.begin();
	for (;it!=vlist.end();it++) {
		delete *it;
	}
	vlist.clear();
}

#else

extern int vftablelist[40];
void addvftable(int a);

// Sicherer ClearPointer mit Typcheckung f�r Elemente die Virtual Function Table haben.
template <class T>
void ClearPointer(std::vector<T *> &vlist) {

	T instanz;
	int vftable_instanz = ((int*) &instanz)[0];

	typename std::vector<T *>::iterator it = vlist.begin();
	for (;it!=vlist.end();it++) {

		int vftable = ((int*) *it)[0];
		if (vftable ==vftable_instanz) {
			delete *it;
			//cout << "ok in " << typeid(instanz).name() << endl;
		} else {
			int i=0;
			while (vftablelist[i]) {
				if (vftablelist[i] == vftable) {
					delete *it;
					break;
				}
				i++;
			}

			if (vftablelist[i]==0) {

				std::cout << "Kaputt in " << typeid(instanz).name() << " " << (int*)vftable << " " <<(int*)vftable_instanz << std::endl;
				//delete *it;

				//	if (typeid(instanz).name()==
			}
		}
	}

	vlist.clear();
}

#endif


/** @ingroup base
 *  Select element number i from pointer vector list.
 *  @param vlist
 *  @param i
 */
template <class T>
T *Select(std::vector<T *> &vlist, int i) {
	int sz = vlist.size();
	return vlist[(i+sz)%sz];
}

/** @ingroup base
 *  Select element number i from pointer vector list.
 *  @param vlist
 *  @param i
 */
template <class T>
T *Select(const std::vector<T *> &vlist, int i) {
	int sz = vlist.size();
	return vlist[(i+sz)%sz];
}

/** @ingroup base
 *  Select element number i from vector list.
 *  @param vlist
 *  @param i
 */
template <class T>
T Select(std::vector<T> &vlist, int i) {
	int sz = vlist.size();
	return vlist[(i+sz)%sz];
}

/** @ingroup base
 * Write vector information to stream.
 * @param os
 * @param vlist
 * @return
 */
template <class T>
std::ostream &operator<<(std::ostream &os, const std::vector<T *> &vlist) {

	typename std::vector<T *>::const_iterator it = vlist.begin();
	for (;it!=vlist.end();it++) {
		os << **it << " ";
	}
	return os;
}

/** @ingroup base
 * Write vector information to stream.
 * @param os
 * @param vlist
 * @return
 */
template <class T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &vlist) {

	typename std::vector<T>::const_iterator it = vlist.begin();
	for (;it!=vlist.end();it++) {
		os << *it << " ";
	}
	return os;
}

/** @ingroup base
 * Write set information to stream.
 * @param os
 * @param vlist
 * @return
 */
template <class T>
std::ostream &operator<<(std::ostream &os, const std::set
	                         <T> &vlist) {

	typename std::set
		<T>::const_iterator it = vlist.begin();
	for (;it!=vlist.end();it++) {
		os << *it << " ";
	}
	return os;
}


/** @ingroup base
 * Write set information to stream.
 * @param os
 * @param vlist
 * @return
 */
template <class T>
std::ostream &operator<<(std::ostream &os, const std::set
	                         <T *> &vlist) {

	typename std::set
		<T *>::const_iterator it = vlist.begin();
	for (;it!=vlist.end();it++) {
		os << **it << std::endl;
	}
	return os;
}



/** @ingroup base
 * Write map information to stream.
 * @param os
 * @param vlist
 * @return
 */
/*template <class T, class U>
std::ostream &operator<<(std::ostream &os, const std::map
	                         <T,U> &vlist) {

	typename std::map
		<T,U>::const_iterator it = vlist.begin();
	for (;it!=vlist.end();it++) {
		os << it->first << " -> " << it->second << endl;
	}
	return os;
}
*/

/** @ingroup base
 * Write map information to stream.
 * @param os
 * @param vlist
 * @return
 */
/*template <class T, class U>
std::ostream &operator<<(std::ostream &os, const std::map
	                         <T,U *> &vlist) {

	typename std::map
		<T,U *>::const_iterator it = vlist.begin();
	for (;it!=vlist.end();it++) {
		os << it->first << " -> " << it->second << endl;
	}
	return os;
}
*/



/** @ingroup base
 * Make deep copy of pointer list.
 * @param source
 * @param dest
 */
template <class T>
void CopyPointer(const std::vector<T *> &source, std::vector<T *> &dest) {

	typename std::vector<T *>::const_iterator it = source.begin();
	for (;it!=source.end();it++) {
		dest.push_back(new T(**it));
	}

}
//@}

}
