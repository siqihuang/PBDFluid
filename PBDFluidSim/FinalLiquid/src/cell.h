#ifndef _CELL_
#define _CELL_

#include "Particle.h"
#include <vector>
class cell
{
public:
	int cell_id;
	int cell_size;
	std::vector<int> cell_particles;
	std::vector<int> cell_neighbors;

	cell(int index, int size);
	~cell();

	void add_particle(int indexP);
	void add_neighbor(int icell);
	void clear_particles();
	std::vector<int> get_particles();
	std::vector<int> get_neighbors();
};
#endif
