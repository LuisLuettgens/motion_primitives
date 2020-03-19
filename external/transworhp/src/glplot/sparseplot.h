#pragma once

#include <GL/glew.h>
#include <GL/gl.h>

#include "baseplot.h"

#include "worhp/C_cs.h"

namespace tw {

class SparsePlot : public BasePlot {
public:
	SparsePlot(WorhpMatrix *d);
	~SparsePlot();
	
protected:

	int rawHighLow(DataStorage &ds) override;
	void drawData(DataStorage &ds, DataStorage &dstop, int ii) override;
	void epsData(DataStorage &ds, DataStorage &dstop, EpsWriter *epsw) override;
	void matlab(std::ostream &os) const override;

	double roundStep(double) const override;
	bool contains(const Point<int> &p) override;
	float MapToY(double d) const override;
	void drawPlotName() const override;

	std::string infostring;
	Point<int> mouse;

private:

	// data worhp matrix
	WorhpMatrix *matrix;
	
	// GL vertices vector
	std::vector<GLfloat> vertices;
	// GL color vector
	std::vector<GLfloat> colors;
	
	// old matrix data from last frame
	std::vector<double> val_old;
	
	GLuint vbo[3];
	
	// vbo initialised?
	bool isInit;
	// nnz worhp matrix last frame
	int nnz_old;
	
	// create VBO
	void init();
	
	//fill matrix data VBO
	void initMatrix();
	
	// clac colors vector
	void updateColor();
};

}
