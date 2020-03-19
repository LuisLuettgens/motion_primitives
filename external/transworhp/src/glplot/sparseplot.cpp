//
// C++ Implementation: plot
//
// Description:
//
//
// Author: Matthias Knauer <knauer@math.uni-bremen.de>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#define GL3_PROTOTYPES 1

#include "sparseplot.h"

#include "xopt_eps.h"

#include "../toolbase/tool.h"

#include <iostream>
#include <iomanip>
#include <algorithm>

enum vboID {matrixVert=0, matrixColor=1, textBox=2};

using namespace std;

namespace tw {

SparsePlot::SparsePlot(WorhpMatrix *d)
        : BasePlot(0), matrix(d),
        vertices(matrix->nnz*4*2,0.1f), colors(matrix->nnz*4*4,1.f/(BasePlot::acsize*.8f*.8f)),
        val_old(matrix->nnz),
        isInit(false),
        nnz_old(0) {

	Low2 = 0.0;
	High2 = static_cast<double>(matrix->nCol); //matrix->nRow;

	Low = 0.0;
	High = static_cast<double>(matrix->nRow); //matrix->nCol;

	Mirror = true;
	
	init();
}

SparsePlot::~SparsePlot() {

	glDeleteBuffers(1,&vbo[matrixVert]);
	glDeleteBuffers(1,&vbo[matrixColor]);
	glDeleteBuffers(1,&vbo[textBox]);
}

void SparsePlot::init() {

	glGenBuffers(1,&vbo[matrixVert]);
	glGenBuffers(1,&vbo[matrixColor]);
	glGenBuffers(1,&vbo[textBox]);
	
	
	GLfloat th = 0.8f;
	
	// vertex test Box
	const vector<GLfloat> box = {
		//box
		0.0f, 0.0f, 0.9f, 			// unten links
		1.0f,236.f/255.f,139.f/255.f,.6f*th,	// Farbe
		1.0f, 0.0f, 0.9f, 			// unten rechts
		1.0f,236.f/255.f,139.f/255.f,.6f*th,	// Farbe
		1.0f,  1.0f, 0.9f, 			// oben rechts
		1.0f,236.f/255.f,139.f/255.f,.9f*th,	// Farbe
		0.0f,  1.0f, 0.9f, 			// oben links
		1.0f,236.f/255.f,139.f/255.f,.9f*th,	// Farbe
		//Rahmen
		0.0f,  0.0f, 1.0f,		// unten links
		0.0f,  0.0f, 0.0f, 0.7f*th, 	// Farbe (schwarz)
		1.0f,  0.0f, 1.0f, 		// unten rechts
		0.0f,  0.0f, 0.0f, 0.7f*th, 	// Farbe (schwarz)
		1.0f,  1.0f, 1.0f, 		// oben rechts
		0.0f,  0.0f, 0.0f, 0.7f*th, 	// Farbe (schwarz)
		0.0f,  1.0f, 1.0f, 		// oben links
		0.0f,  0.0f, 0.0f, 0.7f*th 	// Farbe (schwarz)
	};
	
	glBindBuffer(GL_ARRAY_BUFFER, vbo[textBox]);
	glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*box.size(), box.data(), GL_STATIC_DRAW);
}


void SparsePlot::initMatrix() {
  
	if (!matrix->val) return;
	
	
	if (isInit) {
		glDeleteBuffers(1,&vbo[matrixVert]);
		glGenBuffers(1, &vbo[matrixVert]);
		glBindBuffer(GL_ARRAY_BUFFER, vbo[matrixVert]);
	}
	
  
	vertices.resize(matrix->nnz*4*2,0.1f);
	colors.resize(matrix->nnz*4*4,1.f/(BasePlot::acsize*.8f*.8f));
	
	val_old.resize(matrix->nnz);
	
	//cout << "init nnz:" << matrix->nnz << ", vert.size: " << vertices.size() << " colr.size: " << colors.size() << endl;
	
	for (int i = 0; i < matrix->nnz; i++) {
	  
		int curx, cury, lastx, lasty;

		if (matrix->kind==Matrix_Kind_Vector) {
			curx = 0;
		} else {
			curx = matrix->col[i]-1;
		}
		
		cury = matrix->row[i];


		if (matrix->kind==Matrix_Kind_Vector) {
			lastx = 1;
		} else {
			lastx = matrix->col[i];
		}

		lasty = matrix->row[i]-1;

		//if (lastx-curx<2) lastx = curx+2;
		//if (lasty-cury<2) lasty = cury+2;

		vertices[i*8] = static_cast<GLfloat>(curx);
		vertices[i*8+1] = static_cast<GLfloat>(cury);
		
		vertices[i*8+2] = static_cast<GLfloat>(lastx);
		vertices[i*8+3] = static_cast<GLfloat>(cury);
		
		vertices[i*8+4] = static_cast<GLfloat>(lastx);
		vertices[i*8+5] = static_cast<GLfloat>(lasty);
		
		vertices[i*8+6] = static_cast<GLfloat>(curx);
		vertices[i*8+7] = static_cast<GLfloat>(lasty);
	}
	
	
	glBindBuffer(GL_ARRAY_BUFFER, vbo[matrixVert]);
	glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*vertices.size(), vertices.data(), GL_STATIC_DRAW);
	
	updateColor();
	
	isInit = true;
}

void SparsePlot::updateColor() {
  
	if (isInit) {
		// check if color has changed
		for (size_t i = 0; i < val_old.size(); i++) {
			if (matrix->val[i] != val_old[i]) {
				break;
			} else if (i == val_old.size()-1) {
				return;
			}
		}
	}
  
	//alte werte speichern
	for (size_t i = 0; i < val_old.size(); i++) {
		val_old[i] = matrix->val[i];
	}
  
	for (int i = 0; i < matrix->nnz; i++) {
		
		if (std::abs(matrix->val[i]-1)<1e-6) {
			//SetColor(Green);
			//glColor4f(0.0f,.8f,0.0f,1.0f/sz);
			colors[i*16] = colors[i*16+4] = colors[i*16+8] = colors[i*16+12] = 0.0f;
			colors[i*16+1] = colors[i*16+5] = colors[i*16+9] = colors[i*16+13] = 0.8f;
			colors[i*16+2] = colors[i*16+6] = colors[i*16+10] = colors[i*16+14] = 0.0f;
		}
		else if (std::abs(matrix->val[i]+1)<1e-6) {
			//SetColor(Red);
			//glColor4f(.8f,0.0f,0.0f,1.0f/sz);
			colors[i*16] = colors[i*16+4] = colors[i*16+8] = colors[i*16+12] = 0.8f;
			colors[i*16+1] = colors[i*16+5] = colors[i*16+9] = colors[i*16+13] = 0.0f;
			colors[i*16+2] = colors[i*16+6] = colors[i*16+10] = colors[i*16+14] = 0.0f;
		}
		else if (std::abs(matrix->val[i]+.5)<1e-6) {
			//SetColor(Rose);
			//glColor4f(.8f,.4f,.4f,1.0f/sz);
			colors[i*16] = colors[i*16+4] = colors[i*16+8] = colors[i*16+12] = 0.8f;
			colors[i*16+1] = colors[i*16+5] = colors[i*16+9] = colors[i*16+13] = 0.4f;
			colors[i*16+2] = colors[i*16+6] = colors[i*16+10] = colors[i*16+14] = 0.4f;
		}
		else if (std::abs(matrix->val[i])>1e-6) {
			//SetColor(Cyan);
			//glColor4f(.4f,.4f,1.0f,1.0f/sz);
			colors[i*16] = colors[i*16+4] = colors[i*16+8] = colors[i*16+12] = 0.4f;
			colors[i*16+1] = colors[i*16+5] = colors[i*16+9] = colors[i*16+13] = 0.4f;
			colors[i*16+2] = colors[i*16+6] = colors[i*16+10] = colors[i*16+14] = 1.0f;
		}
		else {
			//SetColor(Grey);
			//glColor4f(.5f,.5f,.5f,1.0f/sz);
			colors[i*16] = colors[i*16+4] = colors[i*16+8] = colors[i*16+12] = 0.5f;
			colors[i*16+1] = colors[i*16+5] = colors[i*16+9] = colors[i*16+13] = 0.5f;
			colors[i*16+2] = colors[i*16+6] = colors[i*16+10] = colors[i*16+14] = 0.5f;
		}
		
	}
	
	//cout << "glBufferData" << endl;
	glBindBuffer(GL_ARRAY_BUFFER, vbo[matrixColor]);
	glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*colors.size(), colors.data(), GL_STATIC_DRAW);
}


double SparsePlot::roundStep(double step) const {
	
	double rd = 1.0;
	double base = 1.0;

	for (int i = 0; i < 10; i++) {

		if (step < base) {
			rd = base;
			break;
		} /*else if (step<2*base) {
			rd = 2*base;
			break;
		} */ else if (step < 5*base) {
			rd = 5*base;
			break;
		}
		base *= 10;
	}
	
	return rd;
}

void SparsePlot::drawData(DataStorage &/*ds*/, DataStorage &/*dstop*/, int /*ii*/) {
	
	if (!isInit || matrix->nnz != nnz_old) {
		// if BFGS change size of HM
		nnz_old = matrix->nnz;
		
		initMatrix();
	}
	
	updateColor();
	
	
	glEnableClientState(GL_COLOR_ARRAY);
	glEnableClientState(GL_VERTEX_ARRAY);
	
	//glColorPointer(4, GL_FLOAT, 0, colors.data());
	//glVertexPointer(2, GL_FLOAT, 0, vertices.data());
	
	glBindBuffer(GL_ARRAY_BUFFER, vbo[matrixVert]);
	glVertexPointer(2, GL_FLOAT, 0, 0);
	
	glBindBuffer(GL_ARRAY_BUFFER, vbo[matrixColor]);
	glColorPointer(4, GL_FLOAT, 0, 0);
	
	glPushMatrix();
	
	if (!IsIcon()) {
		glTranslatef(LBorder, TBorder, 0.0f);
		glScalef((GetWidth()-LBorder-RBorder)/(High2-Low2), -yscale, 1.0f);
		glTranslatef(0.0f, -High, 0.0f);
	} else {
		glScalef((GetWidth())/(High2-Low2), -yscale, 1.0f);
		glTranslatef(0.0f, -High, 0.0f);
	}
	
	// draw
	glDrawArrays(GL_QUADS, 0, matrix->nnz*4);
	
	glPopMatrix();
	
	GLenum err;
	while ((err = glGetError()) != GL_NO_ERROR) {
		cerr << "OpenGL error: " << err << endl;
	}
	
	// deactivate vertex arrays after drawing
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_COLOR_ARRAY);
}

void SparsePlot::drawPlotName() const {
	
	BasePlot::drawPlotName();
	
	if (infostring.size()) {
		
		//int width = GetWidth();
		//int height = GetHeight();

		if (mouseOver > 20 || allplotnames) {
		
			float th = GetHue();

			if (allplotnames) {
				th = 1.0f;
			}

			// Textlaenge
			const float sl = static_cast<float>( Tool::font->StringLength(infostring.c_str()) );

			const float ll = std::min(static_cast<float>(mouse.x+5), GetWidth()-sl-RBorder);
			const float tt = std::max(static_cast<float>(mouse.y-18), static_cast<float>(BBorder));
			
			glEnableClientState(GL_COLOR_ARRAY);
			glEnableClientState(GL_VERTEX_ARRAY);
			
			glBindBuffer(GL_ARRAY_BUFFER, vbo[textBox]);
			glVertexPointer(3, GL_FLOAT, 7*sizeof(GLfloat), 0);
			glColorPointer(4, GL_FLOAT, 7*sizeof(GLfloat), (GLvoid*)(3*sizeof(GLfloat)));
			
			glPushMatrix();
			
			glTranslatef(ll,tt-3,0.0f);
			glScalef(sl+1.0,15.f,1.0f);
			
			// draw
			glDrawArrays(GL_QUADS, 0, 4);
			glDrawArrays(GL_LINE_LOOP, 4, 4);
			
			glPopMatrix();
			
			// deactivate vertex arrays after drawing
			glDisableClientState(GL_VERTEX_ARRAY);
			glDisableClientState(GL_COLOR_ARRAY);
			
			glColor4f(139.f/255.f,129.f/255.f,76.f/255.f,1.0f*th);
			Tool::font->printString(infostring.c_str(),ll+1.0f,tt,1.0f);
		}
	}
}


int SparsePlot::rawHighLow(DataStorage &/*ds*/) {

	Low2 = 0.0;
	High2 = static_cast<double>(matrix->nCol);

	Low = 0.0;
	High = static_cast<double>(matrix->nRow);

	return -1;
}

void SparsePlot::matlab(ostream &os) const {

	if (matrix->kind==Matrix_Kind_LowTri) {

		int lastr = 1, lastc=1;
		
		for (int i=0;i<matrix->nnz;i++) {

			if (matrix->col[i]!=lastc) {

			    while (lastr <= matrix->nRow) {
				os << setw(16) << 0.;
				lastr++;
			    }

			    os << endl;
			    lastc = matrix->col[i];
			    lastr = 1;
			}
			
			while (lastr < matrix->row[i]) {
			    os << setw(16) <<  0.0;
			    lastr++;
			}
			os <<  setw(16) << matrix->val[i] ;
			lastr++;
		}

		while (lastr <= matrix->nRow) {
			os << setw(16) <<  0.0;
			lastr++;
		}

	} else {

		int lastr = 1, lastc=1;
		for (int i=0;i<matrix->nnz;i++) {

			if (matrix->col[i]!=lastc) {

				while (lastr <= matrix->nRow) {
				    os << setw(16) << 0.;
				    lastr++;
				}

				os << endl;
				lastc = matrix->col[i];
				lastr = 1;
			}
			while (lastr < matrix->row[i]) {
				os << setw(16) <<  0.;
				lastr++;
			}
			os <<  setw(16) << matrix->val[i] ;
			lastr++;
		}

		while (lastr <= matrix->nRow) {
			os << setw(16) <<  0.0;
			lastr++;
		}
	}

}


void SparsePlot::epsData(DataStorage &/*ds*/, DataStorage &/*dstop*/, EpsWriter *epsw) {

	for (int i=0;i<matrix->nnz;i++) {

		double curx, cury;
		
		if (matrix->kind==Matrix_Kind_Vector) {
			curx = 0;
		} else {
			curx = epsw->epsMapFloatToX(matrix->col[i]);
		}
		
		cury = epsw->epsMapToY(matrix->row[i]);
		//     double lastx = epsw->epsMapFloatToX(matrix->col[i]);
		//     double lasty = epsw->epsMapToY(matrix->row[i]-1);
		//     if (lastx-curx<1) lastx = curx+1;
		//     if (cury-lasty>1) lasty = cury+1;

		if (fabs(matrix->val[i]-1)<1e-6) {
			epsw->SetLineColorI(23);
		}
		else if (fabs(matrix->val[i]+1)<1e-6) {
			epsw->SetLineColorI(24);
		}
		else if (fabs(matrix->val[i])>1e-6) {
			epsw->SetLineColorI(21);
		}
		else {
			epsw->SetLineColorI(22);
		}

		if (matrix->kind==Matrix_Kind_Vector) {
			epsw->Line((float)epsw->epsMapFloatToX(0.0), (float)cury, (float)epsw->epsMapFloatToX(1), (float)cury);
		}
		else {
			epsw->Dot((float)curx, (float)cury,2);
		}

		// epsw->AddPoint(lastx, lasty);
		// epsw->AddPoint(curx, cury);

		//   epsw->Line();
		// epsw->ClearPoints();
    }
}


float SparsePlot::MapToY(double d) const {

	if (IsIcon()) {
		return static_cast<float>( (High-d)*yscale );
	} else {
		return static_cast<float>( (High-d)*yscale + TBorder );
	}
}



bool SparsePlot::contains(const Point<int> &p)  {

	bool is = p.isin(GetPos(),GetPos()+Point<int>(GetWidth(),GetHeight()));
	
	bool infoflag = false;
	if (is) {
		
		Point<int> tt = p - GetPos();
		
		
		int width = GetWidth();
		float ww = width-LBorder-RBorder;
		float w2 = LBorder;
		int dx = (tt.x - w2 -.5)/ ww *(High2-Low2) +Low2;
		int dy =  High - (tt.y - TBorder)/yscale;
		
		if (dx>=0 && dy>=0 && dx<High2 && dy<High) {
			
			if (matrix->kind==Matrix_Kind_General || matrix->kind==Matrix_Kind_LowTri) {
				for (int i=0;i<matrix->nnz;i++) {
					if (matrix->col[i] == dx+1) {
						if (matrix->row[i] == dy+1) {
							char buf[30];
							sprintf(buf,"[%d/%d] = %.7f",dx,dy, matrix->val[i]);
							infostring = string(buf);
							infoflag = true;
							mouse = tt;
							break;
						}
					}
					//if (matrix->col[i]>dx+1) break;
					
				}
			}
			
			else if (matrix->kind==Matrix_Kind_Vector) {
				for (int i=0;i<matrix->nnz;i++) {
					
						if (matrix->row[i] == dy+1) {
							char buf[30];
							sprintf(buf,"[%d/%d] = %.7f",dx,dy, matrix->val[i]);
							infostring = string(buf);
							infoflag = true;
							mouse = tt;
							break;
						
					}
					//if (matrix->col[i]>dx+1) break;
					
				}
			}
			
		}
		
	}
	
	if (!infoflag) infostring="";
	
	return is;
}

}
