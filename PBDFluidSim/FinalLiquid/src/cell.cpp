#include "cell.h"


cell::cell(int index, int size)
{
	cell_id = index;
	cell_size = size;
}


cell::~cell()
{
}

void cell::add_particle(int indexP){
	cell_particles.push_back(indexP);
}

void cell::add_neighbor(int icell){
	cell_neighbors.push_back(icell);
}

void cell::clear_particles(){
	cell_particles.clear();
}

std::vector<int> cell::get_particles(){
	return cell_particles;
}

std::vector<int> cell::get_neighbors(){
	return cell_neighbors;
}