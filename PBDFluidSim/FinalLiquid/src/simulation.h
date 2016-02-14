#ifndef _SIMULATION_
#define _SIMULATION_

#include "globalHeader.h"
#include "scene.h"
#include "kernel.h"

extern glm::vec3 gravityVector;

class simulation{
private:
public:
	scene *s;
	simulation();
	simulation(scene *s);
	void update(bool GPU_render);
	void updateOnCPU();
	void updateAcc();
	void updateVel();
	void updatePos();
	void collisionDetection();
	void findNeighbor();
	void PBDProjection();
	glm::vec3 spikyKernel(float h,glm::vec3 r);
	float poly6Kernel(float h,glm::vec3 r);
	void calculateLamda(Particle &p);
};

#endif