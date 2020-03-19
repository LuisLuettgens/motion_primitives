/*
 *  GLM library.  Wavefront .obj file format reader/writer/manipulator.
 *
 *  Written by Nate Robins, 1997.
 *  email: ndr@pobox.com
 *  www: http://www.pobox.com/~ndr
 */

#include "model.h"

#include "conversion.h"

#include "../base/vectortools.h"
#include "../base/exception.h"

#include "../core/twstatus.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;

namespace tw {

Model::Model() {}

Model::~Model() {

	vertices.clear();
	normals.clear();
	facetnorms.clear();

	ClearPointer(groups);
	ClearPointer(materials);
}


bool Model::Load(const std::string &path, const std::string &s, double scale, double normalangle) {

	stringstream stat;

	if (s=="")
		return false;
	char buf[100];
	{
		stat << "Loading " << path + s << "... ";
		ifstream is((path+s).c_str());
		pathname = path;

		vector<string> s2 = ToStringArray(s,"/\\");
		s2.pop_back();
		if (s2.size())
			pathname += Join(s2,"/")+ "/";


		if (is) {
			strcpy(buf,(path+s).c_str());
			stat << "ok";

		} else {
			stat << "failed";
			MyStatus("Model", stat.str(), Status::ERR);

			return false;
			/*
						stat.clear();
						stringstream stat2;

						stat2 << "Loading " << s << "... ";
						pathname="";

						vector<string> s2 = ToStringArray(s,"/\\");
						s2.pop_back();
						if (s2.size())
							pathname += Join(s2,"/")+ "/";

						string s3 = s;
						cout << s3.c_str() << "..." <<endl;
						is.open(s3.c_str());
						if (is) {
							strcpy(buf,s3.c_str());
							stat2 << "ok";
						} else {
							stat2 << "failed";

						}
						Status(stat2.str());*/
		}


	}


	glmReadOBJ(buf);

	stat << ", " << groups.size() << " " << "groups";

	glmScale(scale);

	if (normals.size() == 0 || normalangle>0) {
		glmFacetNormals();

		if (normalangle<=0) normalangle=20;

		int r = glmVertexNormals(normalangle);

		stat << ", " << r << " normals generated";
	}

	MyStatus("Model", stat.str(), Status::WARN);

	return true;

}

void Model::Draw(GLuint mode) {

	glPushMatrix();
	glRotatef(90.0f,1.0f,0.0f,0.0f);
//	glmDraw(GLM_SMOOTH|GLM_COLOR);
	glmDraw(mode);
	glPopMatrix();
}

void Model::Draw(GLuint mode, int i) {

	glPushMatrix();
	glRotatef(90.0f,1.0f,0.0f,0.0f);
//	glmDraw(GLM_SMOOTH|GLM_COLOR);
	glmDraw(mode, i);
	glPopMatrix();
}

bool Model::isLoaded() {
	return vertices.size()!=0;
}

size_t Model::countGroups() {
	return groups.size();
}


GLfloat Model::_glmMax(GLfloat a, GLfloat b) {
	if (a > b)
		return a;
	return b;
}


GLfloat Model::_glmAbs(GLfloat f) {
	if (f < 0)
		return -f;
	return f;
}

/*
GLboolean Model::_glmEqual(GLfloat* u, GLfloat* v, GLfloat epsilon) {
	if (_glmAbs(u[0] - v[0]) < epsilon &&
	        _glmAbs(u[1] - v[1]) < epsilon &&
	        _glmAbs(u[2] - v[2]) < epsilon) {
		return GL_TRUE;
	}
	return GL_FALSE;
}


GLfloat* Model::_glmWeldVectors(GLfloat* vectors, GLuint* numvectors, GLfloat epsilon) {
	GLfloat* copies;
	GLuint   copied;
	GLuint   i, j;

	copies = new GLfloat[3*(*numvectors+1)];
	memcpy(copies, vectors, (sizeof(GLfloat) * 3 *(*numvectors + 1)));

	copied = 1;
	for (i = 1; i <= *numvectors; i++) {
		for (j = 1; j <= copied; j++) {
			if (_glmEqual(&vectors[3 * i], &copies[3 * j], epsilon)) {
				goto duplicate;
			}
		}

		/ * must not be any duplicates -- add to the copies array * /
		copies[3 * copied + 0] = vectors[3 * i + 0];
		copies[3 * copied + 1] = vectors[3 * i + 1];
		copies[3 * copied + 2] = vectors[3 * i + 2];
		j = copied;    / * pass this along for below * /
		copied++;

	duplicate:
		/ * set the first component of this vector to point at the correct
		   index into the new copies array * /
		vectors[3 * i + 0] = (GLfloat)j;
	}

	*numvectors = copied-1;
	return copies;
}*/


Model::Group* Model::_glmFindGroup(const std::string &name) {

	std::vector<Group*>::iterator it = groups.begin();
	for (;it!=groups.end();it++) {

		if (name==(*it)->name)
			return (*it);
	}

	return 0;

}


Model::Group* Model::_glmAddGroup(const std::string &name) {
	Group* group;

	group = _glmFindGroup(name);
	if (!group) {
		group = new Group();
		group->name = name;
		groups.push_back(group);
	}

	return group;
}


GLuint Model::_glmFindMaterial(const std::string &name) {
	GLuint i=0;

	//cout << "FindMaterial" << name << endl;

	std::vector<Material*>::iterator it = materials.begin();
	for (;it!=materials.end();it++,i++) {

		if ((*it)->name==name)
			return i;
	}

	/* didn't find the name, so set it as the default material */
	printf("_glmFindMaterial():  can't find material \"%s\".\n", name.c_str());
	i = 0;

	return 0;
}




GLvoid Model::_glmReadMTL(const std::string &name) {
	FILE* file;
	std::string filename;
	char  buf[128];

	filename = pathname + name;


	/* open the file */
	file = fopen(filename.c_str(), "r");
	if (!file) {
		fprintf(stderr, "_glmReadMTL() failed: can't open material file \"%s\".\n",
		        filename.c_str());
		return;
	}


	float a,b,c;
	float alpha = 1;

	/* now, read in the data */
	while (fscanf(file, "%s", buf) != EOF) {
		switch (buf[0]) {
		case '#':    /* comment */
			/* eat up rest of line */
			fgets(buf, sizeof(buf), file);
			break;
		case 'n':    /* newmtl */
			fgets(buf, sizeof(buf), file);
			sscanf(buf, "%s %s", buf, buf);
			materials.push_back(new Material());
			materials.back()->name = buf;
			break;
		case 'N': {
				float a;
				fscanf(file, "%f", &a);

				/* wavefront shininess is from [0, 1000], so scale for OpenGL */
				a *= 128.0/1000.0;

				materials.back()->shininess = a;
			}
			break;
		case 'd':
			fscanf(file, "%f",&a);
			alpha = a;

			break;

		case 'K':
			switch (buf[1]) {
			case 'd':
				fscanf(file, "%f %f %f",&a,&b,&c);
				materials.back()->diffuse = color4(a,b,c,alpha);
				break;
			case 's':
				fscanf(file, "%f %f %f",&a,&b,&c);
				materials.back()->specular = color4(a,b,c,1);

				break;
			case 'a':
				fscanf(file, "%f %f %f",&a,&b,&c);
				materials.back()->ambient = color4(a,b,c,1);

				break;
			default:
				/* eat up rest of line */
				fgets(buf, sizeof(buf), file);
				break;
			}
			break;
		case 'm': /* map_Kd */

			//cout << buf << endl;
			fgets(buf, sizeof(buf), file);


			{
				string imgname(buf);
				imgname = pathname + eatwhitespace(imgname);

		/*		unsigned int a = imgname.find(".jpg");
				if (a!=std::string::npos) {
					imgname.replace(a,4,".png");
				}

				a = imgname.find(".JPG");
				if (a!=std::string::npos) {
					imgname.replace(a,4,".png");
				}*/
				//	cout << " Load Texture for model '" << imgname << "'" << endl;
				materials.back()->texture = new Texture();
				try {
					materials.back()->texture->LoadTexturePNG(imgname,0);
				} catch (Exception &e) {
					//delete materials.back()->texture;
					//materials.back()->texture = 0;
					std::cout << e;
				}

				//	cout << " Loaded" << materials.back().texture->texture << endl;
			}
			break;
		default:
			/* eat up rest of line */
			fgets(buf, sizeof(buf), file);
			break;
		}
	}

	fclose(file);
}



GLvoid Model::_glmPass(FILE* file) {

	Group*    group;      /* current group pointer */
	GLuint    material;   /* current material */
	GLuint    v, n, t;
	char      buf[128];

	group = _glmAddGroup("default");

	float a,b,c;

	material = 0;
	while (fscanf(file, "%s", buf) != EOF) {

		switch (buf[0]) {
			case 'm':
				fgets(buf, sizeof(buf), file);
				sscanf(buf, "%s %s", buf, buf);
				mtllibname = buf;
				_glmReadMTL(buf);
				break;

			case '#':    /* comment */
				/* eat up rest of line */
				fgets(buf, sizeof(buf), file);
				break;
			case 'v':    /* v, vn, vt */
				switch (buf[1]) {
					case '\0':   /* vertex */
						fscanf(file, "%f %f %f",&a,&b,&c);
						vertices.emplace_back(a, b, c);
						break;
					case 'n':    /* normal */
						fscanf(file, "%f %f %f",&a,&b,&c);
						normals.emplace_back(a, b, c);
						break;
					case 't':    /* texcoord */
						fscanf(file, "%f %f",&a,&b);
						//texcoords.emplace_back(a, 1.0f-b);
						texcoords.emplace_back(0.1f, 0.02f);
						break;
				}
				break;
			case 'u':{
				fgets(buf, sizeof(buf), file);
				//cout << "MAT-------" << buf << endl;
				char buf2[100];
				sscanf(buf, " %s", buf2);
				group->material = material = _glmFindMaterial(buf2);
				//cout << "MAT " << buf2 <<endl;
				break;
	}
			case 'g':    /* group */
				/* eat up rest of line */
				fgets(buf, sizeof(buf), file);
				sscanf(buf, "%s", buf);
				group = _glmAddGroup(buf);

	//cout << "GROUP " << buf << endl;
				group->name = buf;
				group->material = material;
				break;
			case 'f':    /* face */
				v = n = t = 0;
				fscanf(file, "%s", buf);
				/* can be one of %d, %d//%d, %d/%d, %d/%d/%d %d//%d */
				if (strstr(buf, "//")) {
					/* v//n */
					sscanf(buf, "%u//%u", &v, &n);

					Triangle a;
					a.v.X() = v;
					a.n.X() = n;
					fscanf(file, "%u//%u", &v, &n);
					a.v.Y() = v;
					a.n.Y() = n;
					fscanf(file, "%u//%u", &v, &n);
					a.v.Z() = v;
					a.n.Z() = n;
					group->triangles.push_back(triangles.size());
					triangles.push_back(a);

					while (fscanf(file, "%u//%u", &v, &n) > 0) {

						Triangle b;
						b.v.X() = a.v.X();
						b.n.X() = a.n.X();
						b.v.Y() = a.v.Z();
						b.n.Y() = a.n.Z();
						b.v.Z() = v;
						b.n.Z() = n;
						group->triangles.push_back(triangles.size());
						triangles.push_back(b);

					}
				} else if (sscanf(buf, "%u/%u/%u", &v, &t, &n) == 3) {
					/* v/t/n */
					Triangle a;
					a.v.X() = v;
					a.t.X() = t;
					a.n.X() = n;
					fscanf(file, "%u/%u/%u", &v, &t, &n);
					a.v.Y() = v;
					a.t.Y() = t;
					a.n.Y() = n;
					fscanf(file, "%u/%u/%u", &v, &t, &n);
					a.v.Z() = v;
					a.t.Z() = t;
					a.n.Z() = n;
					group->triangles.push_back(triangles.size());
					triangles.push_back(a);

					while (fscanf(file, "%u/%u/%u", &v, &t, &n) > 0) {

						Triangle b;

						b.v.X() = a.v.X();
						b.t.X() = a.t.X();
						b.n.X() = a.n.X();
						b.v.Y() = a.v.Z();
						b.t.Y() = a.t.Z();
						b.n.Y() = a.n.Z();
						b.v.Z() = v;
						b.t.Z() = t;
						b.n.Z() = n;
						group->triangles.push_back(triangles.size());
						triangles.push_back(b);
					}

				} else if (sscanf(buf, "%u/%u", &v, &t) == 2) {
					/* v/t */
					Triangle a;
					a.v.X() = v;
					a.t.X() = t;
					fscanf(file, "%u/%u", &v, &t);
					a.v.Y() = v;
					a.t.Y() = t;
					fscanf(file, "%u/%u", &v, &t);
					a.v.Z() = v;
					a.t.Z() = t;
					group->triangles.push_back(triangles.size());
					triangles.push_back(a);

					while (fscanf(file, "%u/%u", &v, &t) > 0) {

						Triangle b;
						b.v.X() = a.v.X();
						b.t.X() = a.t.X();
						b.v.Y() = a.v.Z();
						b.t.Y() = a.t.Z();
						b.v.Z() = v;
						b.t.Z() = t;
						group->triangles.push_back(triangles.size());

						triangles.push_back(b);

					}
				} else {
					/* v */
					Triangle a;

					sscanf(buf, "%u", &v);
					a.v.X() = v;
					fscanf(file, "%u", &v);
					a.v.Y() = v;
					fscanf(file, "%u", &v);
					a.v.Z() = v;
					group->triangles.push_back(triangles.size());
					triangles.push_back(a);


					while (fscanf(file, "%u", &v) > 0) {

						Triangle b;

						b.v.X() = a.v.X();
						b.v.Y() = a.v.Z();
						b.v.Z() = v;
						group->triangles.push_back(triangles.size());
						triangles.push_back(b);


					}
				}
				break;

			default:
				/* eat up rest of line */
				fgets(buf, sizeof(buf), file);
				break;
		}
	}

#if 0
	/* announce the memory requirements */
	printf(" Memory: %d bytes\n",
	       numvertices  * 3*sizeof(GLfloat) +
	       numnormals   * 3*sizeof(GLfloat) *(numnormals ? 1 : 0) +
	       numtexcoords * 3*sizeof(GLfloat) *(numtexcoords ? 1 : 0) +
	       numtriangles * sizeof(GLMtriangle));
#endif

}


GLfloat Model::glmUnitize() {
	GLfloat maxx, minx, maxy, miny, maxz, minz;
	GLfloat cx, cy, cz, w, h, d;
	GLfloat scale;

	/* get the max/mins */
	maxx = minx = vertices[0].X();
	maxy = miny = vertices[0].Y();
	maxz = minz = vertices[0].Z();
	std::vector<Vektor<float> >::iterator it = vertices.begin();
	for (;it!=vertices.end();it++) {

		if (maxx < it->X())
			maxx = it->X();
		if (minx > it->X())
			minx = it->X();

		if (maxy < it->Y())
			maxy = it->Y();
		if (miny > it->Y())
			miny = it->Y();

		if (maxz < it->Z())
			maxz = it->Z();
		if (minz > it->Z())
			minz = it->Z();
	}

	/* calculate model width, height, and depth */
	w = _glmAbs(maxx) + _glmAbs(minx);
	h = _glmAbs(maxy) + _glmAbs(miny);
	d = _glmAbs(maxz) + _glmAbs(minz);

	/* calculate center of the model */
	cx = (maxx + minx) / 2.0;
	cy = (maxy + miny) / 2.0;
	cz = (maxz + minz) / 2.0;

	/* calculate unitizing scale factor */
	scale = 2.0 / _glmMax(_glmMax(w, h), d);

	it = vertices.begin();
	for (;it!=vertices.end();it++) {

		/* translate around center then scale */

		it->X() -= cx;
		it->Y() -= cy;
		it->Z() -= cz;
		it->X() *= scale;
		it->Y() *= scale;
		it->Z() *= scale;
	}

	return scale;
}


GLvoid Model::glmDimensions(GLfloat* dimensions) {
	GLfloat maxx, minx, maxy, miny, maxz, minz;

	/* get the max/mins */
	maxx = minx = vertices[0].X();
	maxy = miny = vertices[0].Y();
	maxz = minz = vertices[0].Z();
	std::vector<Vektor<float> >::iterator it = vertices.begin();
	for (;it!=vertices.end();it++) {

		if (maxx < it->X())
			maxx = it->X();
		if (minx > it->X())
			minx = it->X();

		if (maxy < it->Y())
			maxy = it->Y();
		if (miny > it->Y())
			miny = it->Y();

		if (maxz < it->Z())
			maxz = it->Z();
		if (minz > it->Z())
			minz = it->Z();
	}

	/* calculate model width, height, and depth */
	dimensions[0] = _glmAbs(maxx) + _glmAbs(minx);
	dimensions[1] = _glmAbs(maxy) + _glmAbs(miny);
	dimensions[2] = _glmAbs(maxz) + _glmAbs(minz);
}


GLvoid Model::glmScale(GLfloat scale) {

	std::vector<Vektor<float> >::iterator it = vertices.begin();
	for (;it!=vertices.end();it++) {

		it->X() *= scale;
		it->Y() *= scale;
		it->Z() *= scale;
	}
}


GLvoid Model::glmReverseWinding() {
	GLuint swap;

	std::vector<Triangle>::iterator it = triangles.begin();
	for (;it!=triangles.end();it++) {

		swap = it->v.X();
		it->v.X() = it->v.Z();
		it->v.Z() = swap;

		if (normals.size()) {
			swap = it->n.X();
			it->n.X() = it->n.Z();
			it->n.Z() = swap;
		}

		if (texcoords.size()) {
			swap = it->t.X();
			it->t.X() = it->t.Z();
			it->t.Z() = swap;
		}
	}


	/* reverse facet normals */
	std::vector<Vektor<float> >::iterator it2 = facetnorms.begin();
	for (;it2!=facetnorms.end();it2++) {
		*it2 *= -1;
	}

	/* reverse vertex normals */
	std::vector<Vektor<float> >::iterator it3 = normals.begin();
	for (;it3!=normals.end();it3++) {
		*it3 *= -1;
	}

}


GLvoid Model::glmFacetNormals() {

	/* clobber any old facetnormals */
	facetnorms.clear();

	std::vector<Triangle>::iterator it = triangles.begin();
	for (int i=0;it!=triangles.end();it++,i++) {

		it->findex = i;

		Vektor<float> u,v,w,r;

		u = vertices[ it->v.X() -1];
		v = vertices[ it->v.Y() -1];
		w = vertices[ it->v.Z() -1];

		calcNormal(u,v,w,r);
		facetnorms.push_back(r);

	}
}


int Model::glmVertexNormals(GLfloat angle) {

	struct GLMnode {
		GLuint           index;
		GLboolean        averaged;
		struct GLMnode*  next;
	};

	GLMnode*  node;
	GLMnode*  tail;
	GLMnode** members;

	GLfloat   dot, cos_angle;
	GLuint    i, avg;

	/* calculate the cosine of the angle (in degrees) */
	cos_angle = cos(angle * M_PI / 180.0);

	/* nuke any previous normals */
	normals.clear();

	/* allocate a structure that will hold a linked list of triangle
	   indices for each vertex */
	members = new GLMnode*[vertices.size() + 1];
	for (i = 1; i <= vertices.size(); i++)
		members[i] = NULL;

	/* for every triangle, create a node for each vertex in it */

	std::vector<Triangle>::iterator it = triangles.begin();
	for (int i=0;it!=triangles.end();it++,i++) {

		//	for (i = 0; i < model->numtriangles; i++) {
		node = new GLMnode;
		node->index = i;
		node->next  = members[it->v.X()];
		members[it->v.X()] = node;

		node = new GLMnode;
		node->index = i;
		node->next  = members[it->v.Y()];
		members[it->v.Y()] = node;

		node = new GLMnode;
		node->index = i;
		node->next  = members[it->v.Z()];
		members[it->v.Z()] = node;
	}

	/* calculate the average normal for each vertex */
	for (i = 1; i <= vertices.size(); i++) {
		/* calculate an average normal for this vertex by averaging the
		   facet normal of every triangle this vertex is in */
		node = members[i];
		if (!node)
//			fprintf(stderr, "glmVertexNormals(): vertex w/o a triangle\n");
	MyStatus("Model", "glmVertexNormals(): vertex w/o a triangle", Status::ERR);

		Vektor<float> average(0,0,0);

		avg = 0;
		while (node) {
			/* only average if the dot product of the angle between the two
			   facet normals is greater than the cosine of the threshold
			   angle -- or, said another way, the angle between the two
			   facet normals is less than (or equal to) the threshold angle */

			dot = facetnorms[triangles[node->index].findex]
			      * facetnorms[triangles[members[i]->index].findex];
			if (dot > cos_angle) {
				node->averaged = GL_TRUE;
				average += facetnorms[ triangles[node->index].findex ];
				avg = 1;   /* we averaged at least one normal! */
			} else {
				node->averaged = GL_FALSE;
			}
			node = node->next;
		}

		if (avg) {
			/* normalize the averaged normal */
			average.ReduceToUnit();

			/* add the normal to the vertex normals list */
			normals.push_back(average);
			avg = normals.size();

		}

		/* set the normal of this vertex in each triangle it is in */
		node = members[i];
		while (node) {
			if (node->averaged) {
				/* if this node was averaged, use the average normal */
				if (triangles[node->index].v.X() == i)
					triangles[node->index].n.X() = avg;
				else if (triangles[node->index].v.Y() == i)
					triangles[node->index].n.Y() = avg;
				else if (triangles[node->index].v.Z() == i)
					triangles[node->index].n.Z() = avg;
			} else {
				/* if this node wasn't averaged, use the facet normal */

				normals.push_back(Vektor<float>(facetnorms[ triangles[node->index].findex ]));

				if (triangles[node->index]
				        .v.X() == i)
					triangles[node->index].n.X() = normals.size();
				else if (triangles[node->index]
				         .v.Y() == i)
					triangles[node->index].n.Y() = normals.size();
				else if (triangles[node->index]
				         .v.Z() == i)
					triangles[node->index].n.Z() = normals.size();

			}
			node = node->next;
		}
	}



	/* free the member information */
	for (i = 1; i <= vertices.size(); i++) {
		node = members[i];
		while (node) {
			tail = node;
			node = node->next;
			delete tail;
		}
	}
	delete []members;



	return normals.size();

}


/*
GLvoid Model::glmLinearTexture() {
	GLMgroup *group;
	GLfloat dimensions[3];
	GLfloat x, y, scalefactor;
	GLuint i;

	if (model->texcoords)
		delete []model->texcoords;
	model->numtexcoords = vertices.size();
	model->texcoords=new GLfloat[2*(model->numtexcoords+1)];
	glmDimensions(dimensions);
	scalefactor = 2.0 /
	              _glmAbs(_glmMax(_glmMax(dimensions[0], dimensions[1]), dimensions[2]));

	/ * do the calculations * /
	for (i = 1; i <= vertices.size(); i++) {
		x = vertices[i-1].X() * scalefactor;
		y = vertices[i-1].Z() * scalefactor;
		model->texcoords[2 * i + 0] = (x + 1.0) / 2.0;
		model->texcoords[2 * i + 1] = (y + 1.0) / 2.0;
	}

	/ * go through and put texture coordinate indices in all the triangles * /
	std::vector<GLMgroup*>::iterator it = groups.begin();
	for (;it!=groups.end();it++) {

		for (i = 0; i < (*it)->numtriangles; i++) {
			T((*it)->triangles[i]).tindices[0] = T((*it)->triangles[i]).vindices[0];
			T((*it)->triangles[i]).tindices[1] = T((*it)->triangles[i]).vindices[1];
			T((*it)->triangles[i]).tindices[2] = T((*it)->triangles[i]).vindices[2];
		}

	}

#if 0
	printf("glmLinearTexture(): generated %d linear texture coordinates\n",
	       model->numtexcoords);
#endif
}*/

/*
GLvoid Model::glmSpheremapTexture() {
	GLMgroup* group;
	GLfloat theta, phi, rho, x, y, z, r;
	GLuint i;

	if (model->texcoords)
		delete []model->texcoords;
	model->numtexcoords = normals.size();
	model->texcoords=new GLfloat[2*(model->numtexcoords+1)];

	/ * do the calculations * /
	for (i = 1; i <= normals.size(); i++) {
		z = normals[i-1].X(); / * re-arrange for pole distortion * /
		y = normals[i-1].Y();
		x = normals[i-1].Z();
		r = sqrt((x * x) + (y * y));
		rho = sqrt((r * r) + (z * z));

		if (r == 0.0) {
			theta = 0.0;
			phi = 0.0;
		} else {
			if (z == 0.0)
				phi = M_PI / 2.0;
			else
				phi = acos(z / rho);

#if WE_DONT_NEED_THIS_CODE

			if (x == 0.0)
				theta = M_PI / 2.0; / * asin(y / r); * /
			else
				theta = acos(x / r);
#endif

			if (y == 0.0)
				theta = M_PI / 2.0; / * acos(x / r); * /
			else
				theta = asin(y / r) + (M_PI / 2.0);
		}

		model->texcoords[2 * i + 0] = theta / M_PI;
		model->texcoords[2 * i + 1] = phi / M_PI;
	}

	/ * go through and put texcoord indices in all the triangles * /

	std::vector<GLMgroup*>::iterator it = groups.begin();
	for (;it!=groups.end();it++) {

		for (i = 0; i < (*it)->numtriangles; i++) {
			T((*it)->triangles[i]).tindices[0] = T((*it)->triangles[i]).nindices[0];
			T((*it)->triangles[i]).tindices[1] = T((*it)->triangles[i]).nindices[1];
			T((*it)->triangles[i]).tindices[2] = T((*it)->triangles[i]).nindices[2];
		}

	}

#if 0
	printf("glmSpheremapTexture(): generated %d spheremap texture coordinates\n",
	       model->numtexcoords);
#endif
}
*/
/*
GLvoid Model::glmDelete() {


}*/


Model::Group::~Group() {

	triangles.clear();

}

void Model::glmReadOBJ(const std::string &filename) {

	FILE* file;

	/* open the file */
	file = fopen(filename.c_str(), "r");
	if (!file) {
		fprintf(stderr, "glmReadOBJ() failed: can't open data file \"%s\".\n",
		        filename.c_str());
		return;
	}

#if 0
	/* announce the model name */
	printf("Model: %s\n", filename);
#endif

	_glmPass(file);

	/* close the file */
	fclose(file);

	return;
}


GLvoid Model::glmDraw(GLuint mode) {

	/* do a bit of warning */
	if (mode & GLM_FLAT && !facetnorms.size()) {
		MyStatus("Model", "glmDraw() warning: flat render mode requested "
		                  "with no facet normals defined.",
		         Status::WARN);
		mode &= ~GLM_FLAT;
	}
	if (mode & GLM_SMOOTH && !normals.size()) {
		MyStatus("Model", "glmDraw() warning: smooth render mode requested "
		                  "with no normals defined.",
		         Status::WARN);
		mode &= ~GLM_SMOOTH;
	}
	if (mode & GLM_TEXTURE && !texcoords.size()) {
		MyStatus("Model", "glmDraw() warning: texture render mode requested "
		                  "with no texture coordinates defined.",
		         Status::WARN);
		mode &= ~GLM_TEXTURE;
	}
	if (mode & GLM_FLAT && mode & GLM_SMOOTH) {
		MyStatus("Model", "glmDraw() warning: flat render mode requested "
		                  "and smooth render mode requested (using smooth).",
		         Status::WARN);
		mode &= ~GLM_FLAT;
	}
	if (mode & GLM_COLOR && !materials.size()) {
		MyStatus("Model", "glmDraw() warning: color render mode requested "
		                  "with no materials defined.",
		         Status::WARN);
		mode &= ~GLM_COLOR;
	}
	if (mode & GLM_MATERIAL && !materials.size()) {
		MyStatus("Model", "glmDraw() warning: material render mode requested "
		                  "with no materials defined.",
		         Status::WARN);
		mode &= ~GLM_MATERIAL;
	}
	if (mode & GLM_COLOR && mode & GLM_MATERIAL) {
		MyStatus("Model", "glmDraw() warning: color and material render mode requested "
		                  "using only material mode",
		         Status::WARN);
		mode &= ~GLM_COLOR;
	}
	if (mode & GLM_COLOR)
		glEnable(GL_COLOR_MATERIAL);
	if (mode & GLM_MATERIAL)
		glDisable(GL_COLOR_MATERIAL);

	//glPushMatrix();

	glBegin(GL_TRIANGLES);

	//cout << "GROUPS" << groups.size() << endl;
	std::vector<Group*>::reverse_iterator it = groups.rbegin();
	for (;it!=groups.rend();it++) {

		if (mode & GLM_MATERIAL) {
			glEnd();
			glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT,
			             materials[(*it)->material]->ambient.GetData());
			glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE,
			          materials[(*it)->material]->diffuse.GetData());
			glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR,
			             materials[(*it)->material]->specular.GetData());
			glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS,
			            materials[(*it)->material]->shininess);
			glBegin(GL_TRIANGLES);
		}

		if ((mode & GLM_COLOR)) {// && it!=groups.rbegin()) {
//	cout << materials[(*it)->material]->diffuse ;
//cout << materials[(*it)->material]->diffuse.alpha() << endl;
			glColor4fv(materials[(*it)->material]->diffuse.GetData());
//			glColor3f(1,0,0);

//cout << materials[(*it)->material]->name;
//cout << " " << materials[(*it)->material]->diffuse << endl;

		}

		if (mode & GLM_TEXTURE) {

			if (materials[(*it)->material]->texture) {

				glEnd();

				glEnable(GL_TEXTURE_2D);
				materials[(*it)->material]->texture->Bind();
				//cout << "Binding" << //pathname<<
				//materials[(*it)->material].texture->texture<< endl;

				glBegin(GL_TRIANGLES);

			}
		}

		std::vector<unsigned int>::iterator it2 = (*it)->triangles.begin();
		for (;it2!=(*it)->triangles.end();it2++) {

			Triangle *t = &triangles[ *it2 ];

			if (mode & GLM_FLAT)
				glNormal3fv(facetnorms[ t->findex].data());

			if (mode & GLM_SMOOTH)
				glNormal3fv(normals[t->n.X()-1].data());
			if (mode & GLM_TEXTURE) {
				glTexCoord2fv(texcoords[t->t.X()-1].data());

				//cout << texcoords[t->t.X()-1] << " " ;
				//cout << texcoords[T((*it)->triangles[i]).tindices[0]-1];
			}
			glVertex3fv(vertices[ t->v.X()-1].data());
#if 0

			cout << vertices[T(group->triangles[i]).vindices[0]-1] << endl;
#endif

			if (mode & GLM_SMOOTH)
				glNormal3fv(normals[t->n.Y()-1].data());
			if (mode & GLM_TEXTURE) {
				glTexCoord2fv(texcoords[t->t.Y()-1].data());
				//cout << texcoords[T((*it)->triangles[i]).tindices[1]-1];
			}
			glVertex3fv(vertices[t->v.Y()-1].data());
#if 0

			cout << vertices[T(group->triangles[i]).vindices[1]-1] << endl;
#endif

			if (mode & GLM_SMOOTH)
				glNormal3fv(normals[t->n.Z()-1].data());
			if (mode & GLM_TEXTURE) {
				glTexCoord2fv(texcoords[t->t.Z()-1].data());
				//cout << texcoords[T((*it)->triangles[i]).tindices[2]-1] << endl;
			}

			glVertex3fv(vertices[t->v.Z()-1].data());
#if 0

			cout << vertices[T(group->triangles[i]).vindices[2]-1] << endl;
#endif

		}


		if (mode&GLM_TEXTURE) {

			if (materials[(*it)->material]->texture) {

				glEnd();

				glDisable(GL_TEXTURE_2D);

				glBegin(GL_TRIANGLES);

			}
		}
	}
	glEnd();

	//glPopMatrix();
}




GLvoid Model::glmDraw(GLuint mode, int i) {

	/* do a bit of warning */
	if (mode & GLM_FLAT && !facetnorms.size()) {
		MyStatus("Model", "glmDraw() warning: flat render mode requested "
		                  "with no facet normals defined.",
		         Status::WARN);
		mode &= ~GLM_FLAT;
	}
	if (mode & GLM_SMOOTH && !normals.size()) {
		MyStatus("Model", "glmDraw() warning: smooth render mode requested "
		                  "with no normals defined.",
		         Status::WARN);
		mode &= ~GLM_SMOOTH;
	}
	if (mode & GLM_TEXTURE && !texcoords.size()) {
		MyStatus("Model", "glmDraw() warning: texture render mode requested "
		                  "with no texture coordinates defined.",
		         Status::WARN);
		mode &= ~GLM_TEXTURE;
	}
	if (mode & GLM_FLAT && mode & GLM_SMOOTH) {
		MyStatus("Model", "glmDraw() warning: flat render mode requested "
		                  "and smooth render mode requested (using smooth).",
		         Status::WARN);
		mode &= ~GLM_FLAT;
	}
	if (mode & GLM_COLOR && !materials.size()) {
		MyStatus("Model", "glmDraw() warning: color render mode requested "
		                  "with no materials defined.",
		         Status::WARN);
		mode &= ~GLM_COLOR;
	}
	if (mode & GLM_MATERIAL && !materials.size()) {
		MyStatus("Model", "glmDraw() warning: material render mode requested "
		                  "with no materials defined.",
		         Status::WARN);
		mode &= ~GLM_MATERIAL;
	}
	if (mode & GLM_COLOR && mode & GLM_MATERIAL) {
		MyStatus("Model", "glmDraw() warning: color and material render mode requested "
		                  "using only material mode",
		         Status::WARN);
		mode &= ~GLM_COLOR;
	}
	if (mode & GLM_COLOR)
		glEnable(GL_COLOR_MATERIAL);
	if (mode & GLM_MATERIAL)
		glDisable(GL_COLOR_MATERIAL);

	//glPushMatrix();

	glBegin(GL_TRIANGLES);

	//cout << "GROUPS" << groups.size() << endl;
	std::vector<Group*>::iterator it = groups.begin();
	it += i;
	//for (;it!=groups.rend();it++)
	{

		if (mode & GLM_MATERIAL) {
			glEnd();
			glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT,
			             materials[(*it)->material]->ambient.GetData());
			glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE,
			          materials[(*it)->material]->diffuse.GetData());
			glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR,
			             materials[(*it)->material]->specular.GetData());
			glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS,
			            materials[(*it)->material]->shininess);
			glBegin(GL_TRIANGLES);
		}

		if ((mode & GLM_COLOR)) {// && it!=groups.rbegin()) {
//	cout << materials[(*it)->material]->diffuse ;
//cout << materials[(*it)->material]->diffuse.alpha() << endl;
			glColor4fv(materials[(*it)->material]->diffuse.GetData());
//			glColor3f(1,0,0);

//cout << materials[(*it)->material]->name;
//cout << " " << materials[(*it)->material]->diffuse << endl;

		}

		if (mode & GLM_TEXTURE) {

			if (materials[(*it)->material]->texture) {

				glEnd();

				glEnable(GL_TEXTURE_2D);
				materials[(*it)->material]->texture->Bind();
				//cout << "Binding" << //pathname<<
				//materials[(*it)->material].texture->texture<< endl;

				glBegin(GL_TRIANGLES);

			}
		}

		std::vector<unsigned int>::iterator it2 = (*it)->triangles.begin();
		for (;it2!=(*it)->triangles.end();it2++) {

			Triangle *t = &triangles[ *it2 ];
			
			if (mode & GLM_FLAT)
				glNormal3fv(facetnorms[ t->findex].data());

			if (mode & GLM_SMOOTH)
				glNormal3fv(normals[t->n.X()-1].data());
			if (mode & GLM_TEXTURE) {
				glTexCoord2fv(texcoords[t->t.X()-1].data());

				//cout << texcoords[t->t.X()-1] << " " ;
				//cout << texcoords[T((*it)->triangles[i]).tindices[0]-1];
			}
			glVertex3fv(vertices[ t->v.X()-1].data());
#if 0

			cout << vertices[T(group->triangles[i]).vindices[0]-1] << endl;
#endif

			if (mode & GLM_SMOOTH)
				glNormal3fv(normals[t->n.Y()-1].data());
			if (mode & GLM_TEXTURE) {
				glTexCoord2fv(texcoords[t->t.Y()-1].data());
				//cout << texcoords[T((*it)->triangles[i]).tindices[1]-1];
			}
			glVertex3fv(vertices[t->v.Y()-1].data());
#if 0

			cout << vertices[T(group->triangles[i]).vindices[1]-1] << endl;
#endif

			if (mode & GLM_SMOOTH)
				glNormal3fv(normals[t->n.Z()-1].data());
			if (mode & GLM_TEXTURE) {
				glTexCoord2fv(texcoords[t->t.Z()-1].data());
				//cout << texcoords[T((*it)->triangles[i]).tindices[2]-1] << endl;
			}

			glVertex3fv(vertices[t->v.Z()-1].data());
#if 0

			cout << vertices[T(group->triangles[i]).vindices[2]-1] << endl;
#endif

		}


		if (mode&GLM_TEXTURE) {

			if (materials[(*it)->material]->texture) {

				glEnd();

				glDisable(GL_TEXTURE_2D);

				glBegin(GL_TRIANGLES);

			}
		}
	}
	glEnd();

	//glPopMatrix();
}





/*GLuint Model::glmList(GLuint mode) {
	GLuint list;

	list = glGenLists(1);
	glNewList(list, GL_COMPILE);
	glmDraw(mode);
	glEndList();

	return list;
}*/

#if 0
GLvoid Model::glmWeld(GLfloat epsilon) {

	GLfloat* vectors;
	GLfloat* copies;
	GLuint   numvectors;
	GLuint   i;

	/* vertices */
	numvectors = model->numvertices;
	vectors    = model->vertices;
	copies = _glmWeldVectors(vectors, &numvectors, epsilon);

	printf("glmWeld(): %d redundant vertices.\n",
	       model->numvertices - numvectors - 1);

	for (i = 0; i < model->numtriangles; i++) {
		T(i).vindices[0] = (GLuint)vectors[3 * T(i).vindices[0] + 0];
		T(i).vindices[1] = (GLuint)vectors[3 * T(i).vindices[1] + 0];
		T(i).vindices[2] = (GLuint)vectors[3 * T(i).vindices[2] + 0];
	}

	/* free space for old vertices */
	delete []vectors;

	/* allocate space for the new vertices */
	model->numvertices = numvectors;
	model->vertices = new GLfloat[3 *(model->numvertices + 1)];

	/* copy the optimized vertices into the actual vertex list */
	for (i = 1; i <= model->numvertices; i++) {
		model->vertices[3 * i + 0] = copies[3 * i + 0];
		model->vertices[3 * i + 1] = copies[3 * i + 1];
		model->vertices[3 * i + 2] = copies[3 * i + 2];
	}

	delete []copies;


}

#endif


Model::Material::Material() {

	name = "";
	shininess = 0;
	diffuse = color4(0.8,0.8,0.8,1.0);
	ambient = color4(0.2,0.2,0.2,1.0);
	specular = color4(0.0,0.0,0.0,1.0);

	texture = 0;

}

Model::Material::~Material() {

	delete texture;

}

}
