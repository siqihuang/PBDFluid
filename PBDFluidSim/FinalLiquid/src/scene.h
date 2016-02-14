#ifndef _SCENE_
#define _SCENE_

#include "globalHeader.h"
#include "container.h"
#include "grid.h"
#include "io.h"
#include "primitive.h"

class scene{
private:
public:
	vector<Particle> particles; 
	container *cont;
	vector<primitive*> pri;

	scene();
	bool init(fileInput *fileinput,container *cont);
	void drawPrimitive(const VBO& vbos);
	void drawParticle(const VBO& vbos,bool GPU_render);
};

#endif