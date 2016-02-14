#ifndef _CONTAINER_
#define _CONTAINER_

#include "globalHeader.h"
#include "Particle.h"
#include "grid.h"

class container{
protected:
	glm::vec3 pos;//containner position
	glm::vec3 boundingBox;//container boundingBox
	float cell_size;//containner grid cell size
public:
	glm::vec3 minPos;
	grid *g;
	container();
	container(glm::vec3 pos,float cell_size);
	virtual void staticIntersection(Particle &p);
	virtual void draw();
};

class container_cube:public container{
private:
	glm::vec3 half_size;
	glm::vec3 lower_bound,upper_bound;
public:
	container_cube(glm::vec3 pos,glm::vec3 dimension,float cell_size);
	virtual void staticIntersection(Particle &p);
	virtual void draw();
};

class container_sphere:public container{
private:
	float radius;
public:
	container_sphere(glm::vec3 pos,float radius,float cell_size);
	virtual void staticIntersection(Particle &p);
	virtual void draw();
};

#endif