#ifndef _KERNEL_
#define _KERNEL_

#include <cuda.h>
#include <cuda_runtime.h>
#include "globalHeader.h"
#include "kdtree.h"

extern glm::vec3 gravityVector;

struct GPUParticle{
	glm::vec3 pos,vel,acc,newPos;//position, velocity, acceleration and new position
	glm::vec3 cellId;
	float lamda;
	float mass;
	int neighbor[450];//this num must be consistent with GPU_particle_neighbor_limit
	int index;//index of current recorded neighbor
};

struct GPUCell{
	glm::vec3 pos;
	int particles[50];//this num must be consistent with GPU_cell_particle_limit
	int index;//index of current recorded particle
};

struct GPUContainer{
	int container_type;
	glm::vec3 cubeDim;
	float sphereRadius;
	glm::vec3 minPos;
	glm::vec3 boundingBox;
	glm::vec3 cell_dim;
	float cell_size;
	GPUCell *cell;
};

struct GPUPrimitive{
	int primitive_type;
	glm::vec3 pos;
	glm::vec3 cubeDim;
	float sphereRadius;
	kdtree *mesh;
	glm::vec3 *objVertex;
	glm::vec3 *objNormal;
	int *objIndices;
};

kdtree *initTree(kdtree *root);
void initParticleOnGPU(glm::vec3 particleDim,glm::vec3 particleMinPos,float particleDis);
void initCubeContainerOnGPU(glm::vec3 dimension,float cell_size);
void initSphereContainerOnGPU(float radius,float cell_size);
void initCellOnGPU(glm::vec3 boundingBox,glm::vec3 cell_dim,float cell_size);
void initPrimitiveOnGPU(GPUPrimitive *primitive,int num);
void updateOnGPU();
void updateAccOnGPU();
void updateVelOnGPU();
void updatePosOnGPU();
void updateStateOnGPU();
void putParticleInCellOnGPU();
void clearCellParticleOnGPU();
void clearParticleNeighborOnGPU();
void containerCollisionOnGPU();
void primitiveCollisionOnGPU();
void findNeighborOnGPU();
void PBDProjectionOnGPU();
void calculateLamdaOnGPU();
void calculateNewVelOnGPU();
glm::vec3 spikyKernelOnGPU(float h,glm::vec3 r);
float poly6KernelOnGPU(float h,glm::vec3 r);
glm::vec3 getPos(int n);
void deleteData();

#endif