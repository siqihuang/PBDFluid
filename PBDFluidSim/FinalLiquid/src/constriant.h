#ifndef _CONSTRIANT_
#define _CONSTRIANT_

const float PLANE_COLOR[3] = {0.2, 0.2, 0.2}; //grey
const float dt = 1.0 / 100.0;
const float gravity = -980;
const float rebound_rest = 0.6;
const float particle_radius=0.6;
const float PI=3.1415926;
const float smooth_radius = 1.5;
const float rest_density = 300.0f;
const float relaxation = 0.1f;
const float Epsilon = 1e-3;
const float colision_epsilon = 1e-3;
const float GPU_particle_neighbor_limit=450;
const float GPU_cell_particle_limit=50;

#endif