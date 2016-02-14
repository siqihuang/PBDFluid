#ifndef _PRIMITIVE_
#define _PRIMITIVE_

#include "globalHeader.h"
#include "openglHeader.h"
#include "Particle.h"
#include "kernel.h"
#include <fstream>

class primitive{
protected:
	glm::vec3 pos;
	int primitive_type;
	vector<glm::vec3> m_positions,m_normals,m_colors;
	vector<unsigned short> m_indices;

public:
	primitive();
	primitive(glm::vec3 pos);
	void draw(const VBO& vbos);
	virtual void init();
	virtual void staticIntersectionTest(Particle &p);
	virtual GPUPrimitive getGPUPrimitive();
	kdtree tree;
};

class cube:public primitive{
private:
	glm::vec3 dimension,m_hf_dims,lower_bound,upper_bound;
public:
	cube();
	cube(glm::vec3 pos,glm::vec3 dimension);
	virtual void init();
	virtual void staticIntersectionTest(Particle &p);
	virtual GPUPrimitive getGPUPrimitive();
};

class sphere:public primitive{
private:
	float radius;
public:
	sphere();
	sphere(glm::vec3 pos,float radius);
	virtual void init();
	virtual void staticIntersectionTest(Particle &p);
	virtual GPUPrimitive getGPUPrimitive();
};

class objmesh:public primitive{
private:
	float m_scaling;
public:
	objmesh();
	objmesh(glm::vec3 pos,float m_scaling);
	virtual void init();
	virtual void staticIntersectionTest(Particle &p);
	void objmesh::read_from_file(char* filename);
	virtual GPUPrimitive getGPUPrimitive();
};

#endif