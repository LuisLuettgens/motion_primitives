//
// Author: Matthias Knauer <knauer@math.uni-bremen.de>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
#ifndef dynloader_h
#define dynloader_h

#include <string>
#include <map>

/** @ingroup base
 * Function handle to create a new instance of a class.
 * 
 * Used by class DynLoader.
 */
typedef void* (*newinstance) ();

/** @ingroup base
 *  Macro to add a class to the registration system. 
 */
#define REGISTER_CLASS(A) \
static void* newClass() { \
	return new A();\
} \
static void registerClass() { \
	dynloader.registerClass(#A, A::newClass);\
}

/** @ingroup base
 *  Macro to add a OpenGL class to the registration system. 
 */
#define REGISTER_GL_CLASS(A,B) \
static void* newClass() { \
	return new A();\
} \
static void registerClass() { \
	dynloader.registerClass(B, A::newClass);\
}

// Anm.: #A erstellt aus dem übergegebenen Parameter A den String "A".

/** @ingroup base
 *  @brief Create instances of classes by class name.
 *
 * Add the REGISTER_CLASS()-Macro to the class declaration.
 * Register a class by calling the static function Classname::registerClass(); at program start.
 * dynloader::newClass("Classname") then returns a new instance of this class.
 */
class DynLoader {
public:

	/**
	 * Constructor.
	 */
	DynLoader() {}
	/**
	 * Destructor.
	 */
	~DynLoader();
	/**
	 * Register a class.
	 * @param name Name of the class.
	 * @param a Handle to function returning a new instance of the class.
	 */
	void registerClass(const std::string &name, newinstance a);
	/**
	 * Create instance of class.
	 * @param name Name of the class.
	 * @return New instance.
	 */
	void *newClass(const std::string &name);

	//friend std::ostream &operator<<(std::ostream &os, const DynLoader &s);

private:
	//std::map<std::string, newinstance, std::less<std::string>, std::__malloc_alloc_template<0> > m;
	std::map<std::string, newinstance> m;

};

/** Global instance of DynLoader. */
extern DynLoader dynloader;

#endif
