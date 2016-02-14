#ifndef _GRID_
#define _GRID_

#include "globalHeader.h"
#include "cell.h"
using namespace std;

class grid
{
public:
	float grid_length;
	float grid_height;
	float grid_depth;
	float cell_size;
	glm::vec3 gdim;

	grid(glm::vec3 boundingBox, float cell_size);
	~grid();
	
	vector<vector<vector<cell>>> cells;

	vector<int> neighbor_particles(int particleID, glm::vec3 cellID);

	float get_length();
	float get_height();
	float get_depth();
	void putParticle(Particle &p,glm::vec3 offset,int particleId);
	void clearCellParticle();
	void print();
};
#endif