#pragma once
#include <GL/glew.h>
#include <GL/gl.h>

#include "model.h"

#include "xmlio.h"

#include "../base/vektor.h"

#include <string>

#ifndef DllExport
#ifdef _WIN32
#define DllExport __declspec(dllexport)
#else
#define DllExport
#endif
#endif

#ifdef _MSC_VER
#pragma warning(disable : 4251)
#endif

namespace tw {

class DllExport glObject {
public:
	glObject();
	~glObject();

	void Init(XMLNode *n);
	void InitObject(const std::string &path);
	void InitSplitObject(const std::string &path);
	void Draw(Vektor<float>&pos, float phi, float psi, int i=1);

	int countObj() {return splitObjects;}

private:

	GLuint list;
	int splitObjects;
	double scale;
	std::string filename;
	std::string id;
	Model *m2;
	double normalangle;
};

}
