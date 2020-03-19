#include "dynloader.h"

#include <iostream>

#include "exception.h"

using namespace std;

DynLoader dynloader;

DynLoader::~DynLoader() {
	m.clear();
}


void *DynLoader::newClass(const std::string &name) {

	newinstance a = m[name];

	if (a) {
		//cout << "NEW CLASS " << name << endl;
		return a();
	} else {
		Exception e;
		e << "Class " << name << " not found!";
		throw e;
	}
	return 0;
}


void DynLoader::registerClass(const std::string &name, newinstance a) {

	m[name] = a;
	
}

/*
std::ostream &operator<<(std::ostream &os, const DynLoader &s) {
	os << s.m;
	return os;
}
*/
