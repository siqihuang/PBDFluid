#ifndef _PARTICLE_
#define _PARTICLE_

#include "globalHeader.h"
#include "openglHeader.h"

class Particle{
private:
	float radius,rest;
public:
	glm::vec3 pos,vel,acc,newPos;//position, velocity, acceleration and new position
	glm::vec3 dist;//dist inside primitive
	glm::vec3 cellId;

	/*
	//used for shader to draw
	*/
	std::vector<glm::vec3> m_positions,m_colors,m_normals;
	std::vector<unsigned short> m_indices;
	/*
	//used for shader to draw
	*/

	float density;
	float lamda;
	float mass;
	std::vector<int> neighbor;

	Particle();
	Particle(float radius);
	Particle(float radius,glm::vec3 pos);
	void draw(const VBO& vbos);
	void init();
	float getRadius();
	float getRest();
};

#endif