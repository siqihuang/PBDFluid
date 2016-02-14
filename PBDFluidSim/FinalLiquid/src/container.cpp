#include "container.h"

container::container(){}

container::container(glm::vec3 pos,float cell_size){
	this->pos=pos;
	this->cell_size=cell_size;
}

void container::staticIntersection(Particle &p){
	
}

void container::draw(){}

container_cube::container_cube(glm::vec3 pos,glm::vec3 dimension,float cell_size):container(pos,cell_size){
	this->half_size=dimension/2.0f;
	boundingBox=dimension;

	lower_bound=pos-half_size;
	upper_bound=pos+half_size;
	this->minPos=glm::vec3(0,0,0)-boundingBox/2.0f;//need improvement

	g=new grid(boundingBox,cell_size);
}

void container_cube::staticIntersection(Particle &p){
	glm::vec3 pos=p.newPos;
	float radius=p.getRadius();
	for(int i=0;i<3;++i){
		if(pos[i]-radius-lower_bound[i]<=0){
			p.newPos[i]+=(radius+lower_bound[i]-pos[i]);
			p.newPos[i]+=(0.04*1.0*rand()/(RAND_MAX+1)-0.02);
			p.vel[i]*=-rebound_rest;
			p.newPos[i]+=p.vel[i]*dt/2;
		}
		if(pos[i]+radius-upper_bound[i]>=0){
			p.newPos[i]-=(pos[i]+radius-upper_bound[i]);
			p.newPos[i]+=(0.04*1.0*rand()/(RAND_MAX+1)-0.02);
			p.vel[i]*=-rebound_rest;
			p.newPos[i]+=p.vel[i]*dt/2;
		}
	}
}

void container_cube::draw(){
	glColor3fv(PLANE_COLOR);
	glScaled(half_size[0]*2,half_size[1]*2,half_size[2]*2);
	glutWireCube(1.0);
	glPopMatrix();
}

container_sphere::container_sphere(glm::vec3 pos,float radius,float cell_size):container(pos,cell_size){
	this->radius=radius;
	this->minPos=glm::vec3(0,0,0)-glm::vec3(radius,radius,radius);//need improvement
	boundingBox=2.0f*glm::vec3(radius);

	g=new grid(boundingBox,cell_size);
}

void container_sphere::staticIntersection(Particle &p){
	glm::vec3 pos=p.newPos;
	float radius=p.getRadius();
	float dis=glm::length(pos-this->pos);
	if(dis+radius>this->radius){
		dis=dis+radius-this->radius;
		glm::vec3 nor=glm::normalize(pos-this->pos);
		p.newPos-=nor*dis;

		p.newPos+=(float)(0.04*1.0*rand()/(RAND_MAX+1)-0.02)*nor;
		p.newPos-=glm::dot(p.vel,nor)*nor*dt/2.0f;

		glm::vec3 tmp=p.vel-glm::dot(p.vel,nor)*nor;
		p.vel=-glm::dot(p.vel,nor)*nor*rebound_rest+rebound_rest*tmp;
	}
}

void container_sphere::draw(){
	glColor3fv(PLANE_COLOR);
	glScaled(radius,radius,radius);
	glutWireSphere(1,20,20);
	glPopMatrix();
}