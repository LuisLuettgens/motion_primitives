#include "globject.h"

namespace tw {

glObject::glObject() : list(0), scale(1.0), m2(nullptr), normalangle(0.0) {}


void glObject::Init(XMLNode *n) {

	if (n) {
		filename = n->GetText();
		std::string s = n->GetAttribute("scale");
		if (!s.empty()) {
			scale = std::stod(s);
		}

		id = n->GetAttribute("id");

		s = n->GetAttribute("normal");
		if (!s.empty()) {
			normalangle = std::stod(s);
		}
	}
}


glObject::~glObject() {

	glDeleteLists(list,2);

	delete m2;
}


void glObject::Draw(Vektor<float>&pos, float phi, float psi, int i) {

	glPushMatrix();
	glTranslatef(pos.X(),pos.Y(),pos.Z());
	glRotatef(phi,0,0,1);
	glRotatef(psi,0,1,0);
	glCallList(list + i);
	//m2->Draw(GLM_SMOOTH|GLM_MATERIAL|GLM_TEXTURE,i); // GLM_MATERIAL
		
	glPopMatrix();
}


void glObject::InitSplitObject(const std::string &path) {

/*
	if (m.isLoaded()) {
		glNewList(list + 1 , GL_COMPILE);
		m.Draw();
		glEndList();
	}
	else if (m.Load(path,filename,scale)) {
		glNewList(list + 1 , GL_COMPILE);
		m.Draw();
		glEndList();
	}
	*/
	
		//Model m2;
	m2 = new Model();
	splitObjects = 0;
	
	if (m2->Load(path + "modell/",filename,scale,normalangle)) {

		splitObjects = m2->countGroups();
		list = glGenLists(splitObjects);

		for (int i = 0; i < splitObjects; i++) {
			glNewList(list + i , GL_COMPILE);
			m2->Draw(GLM_SMOOTH|GLM_MATERIAL|GLM_TEXTURE,i); // GLM_MATERIAL
			glEndList();
		}
	}
}

void glObject::InitObject(const std::string &path) {

	list = glGenLists(3);
/*
	if (m.isLoaded()) {
		glNewList(list + 1 , GL_COMPILE);
		m.Draw();
		glEndList();
	}
	else if (m.Load(path,filename,scale)) {
		glNewList(list + 1 , GL_COMPILE);
		m.Draw();
		glEndList();
	}
	*/
	
		//Model m2;
	m2 = new Model();
	
	if (m2->Load(path + "modell/",filename,scale,normalangle)) {
		
		glNewList(list + 0, GL_COMPILE);
		m2->Draw(GLM_SMOOTH|GLM_MATERIAL|GLM_TEXTURE); // GLM_MATERIAL
		glEndList();

		glNewList(list + 1, GL_COMPILE);
		m2->Draw(GLM_SMOOTH|GLM_COLOR|GLM_TEXTURE); // GLM_MATERIAL
		glEndList();
	
		glNewList(list + 2, GL_COMPILE);
		m2->Draw(GLM_SMOOTH  /*|GLM_COLOR*/);
		glEndList();
	}
}

}
