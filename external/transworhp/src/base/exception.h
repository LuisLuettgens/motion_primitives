//
// Author: Matthias Knauer <knauer@math.uni-bremen.de>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
#ifndef exception_h
#define exception_h

#include "../base/defines.h"

#include <string>
#include <sstream>

/** @defgroup base Basic Structures
 * @brief Collection of general classes for programming.
 */

/** @ingroup base
 *  @brief Streamable exception class.
 *
 *  Stream text into the Exception and throw it, or just make a warning.
 */
class Exception {

public:
	/** Constructor. */
	Exception() {}
	/** Constructor. */
	Exception(const std::string &e) {
		what << e;
	}
	/** Copy-Constructor. */
	Exception(const Exception &e) {
		what << e.what.str();
	}

	/** Print out Exception information. */
	friend std::ostream &operator<<(std::ostream &os, const Exception &e);

	/** Anything /const/ can be added to Exception stream. */
	template <class T>
	friend Exception &operator<<(Exception &e, const T &d) {
		e.what<<d;
		return e;
	}
	
	/**
	 * Write Exception as a warning to cout.
	 * @param col Color of warning text.
	 *     - 0: black
	 *     - 1: red
	 *     - 2: green
	 *     - 3: blue
	 *     - 4: grey
	 *     - 5: dark red
	 */
	void warn(int col=3);


private:
	/** Exception information. */
	std::stringstream what;
};

#endif
