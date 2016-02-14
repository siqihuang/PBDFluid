#ifndef _IOINPUT
#define _IOINPUT

#include "globalHeader.h"
#include "openglHeader.h"
#include "primitive.h"
#include <fstream>

class fileInput{
public:
	fileInput();
	fileInput(const char* filename);
	bool readFile();
	glm::vec3 getCubeDim();
	float getSphereRadius();
	glm::vec3 getParticleDim();
	glm::vec3 getParticleMinPos();
	glvuVec3f getCameraLookAt();
	glvuVec3f getCameraPos();
	glvuVec3f getCameraUp();
	float getParticleDis();
	int getContainerType();
	bool isGPURender();

	const char* filename;
	std::ifstream input;
	vector<primitive*> pri;

private:
	int container_type;
	float sphere_radius;
	glm::vec3 cube_dim;
	glm::vec3 particle_dim,particle_minPos;
	glvuVec3f camera_lookat,camera_pos,camera_up;
	float particle_dis;
	bool GPU_render;
};

#endif