#include "scene.h"
#include "kernel.h"

scene::scene(){}

bool scene::init(fileInput *fileinput,container *cont){
	this->cont=cont;
	this->pri=fileinput->pri;

	glm::vec3 particle_dim=fileinput->getParticleDim();
	glm::vec3 particle_minPos=fileinput->getParticleMinPos();
	float dis=fileinput->getParticleDis();

	for(int i=0;i<particle_dim.x;++i){
		for(int j=0;j<particle_dim.y;++j){
			for(int k=0;k<particle_dim.z;++k){
				glm::vec3 pos=particle_minPos+glm::vec3(i,j,k)*dis;
				Particle p(particle_radius,pos);
				particles.push_back(p);
			}//k
		}//j
	}//i

	return true;
}

void scene::drawPrimitive(const VBO& vbos){
	for(int i=0;i<pri.size();++i){
		pri[i]->draw(vbos);
	}
}

void scene::drawParticle(const VBO& vbos,bool GPU_render){
	if(GPU_render){
		for(int i=0;i<particles.size();++i){
			glm::vec3 pos=getPos(i);
			glPushMatrix();
			glTranslated(pos[0],pos[1],pos[2]);
			glutSolidSphere(particle_radius, 10, 10);
			glPopMatrix();
			//cout<<"outsize "<<pos.x<<","<<pos.y<<","<<pos.z<<endl;
		}
	}
	else{
		for(int i=0;i<particles.size();++i){
			particles[i].draw(vbos);
		}
	}
}