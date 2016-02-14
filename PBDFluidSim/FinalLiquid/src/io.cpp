#include "io.h"
fileInput::fileInput(){}

fileInput::fileInput(const char *filename){
	this->filename=filename;
}

bool fileInput::readFile(){
	input.open(filename);
	std::string s;
	while(!input.eof()){
		input>>s;
		if(s=="CONTAINER_CUBE"){
			container_type=CONTAINER::CONTAINER_CUBE;
			input>>cube_dim.x;
			input>>cube_dim.y;
			input>>cube_dim.z;
		}
		else if(s=="CONTAINER_SPHERE"){
			container_type=CONTAINER::CONTAINER_SPHERE;
			input>>sphere_radius;
		}
		else if(s=="EYE"){
			input>>camera_pos.x;
			input>>camera_pos.y;
			input>>camera_pos.z;
		}
		else if(s=="LOOKAT"){
			input>>camera_lookat.x;
			input>>camera_lookat.y;
			input>>camera_lookat.z;
		}
		else if(s=="UP"){
			input>>camera_up.x;
			input>>camera_up.y;
			input>>camera_up.z;
		}
		else if(s=="PARTICLE"){
			input>>particle_dim.x;
			input>>particle_dim.y;
			input>>particle_dim.z;
		}
		else if(s=="MINPOS"){
			input>>particle_minPos.x;
			input>>particle_minPos.y;
			input>>particle_minPos.z;
		}
		else if(s=="DIS"){
			input>>particle_dis;
		}
		else if(s=="CUBE"){
			input>>s;
			glm::vec3 pos,dim;
			input>>pos.x;input>>pos.y;input>>pos.z;
			input>>s;
			input>>dim.x;input>>dim.y;input>>dim.z;
			primitive *p=new cube(pos,dim);
			pri.push_back(p);
		}
		else if(s=="SPHERE"){
			input>>s;
			glm::vec3 pos;
			float radius;
			input>>pos.x;input>>pos.y;input>>pos.z;
			input>>s;
			input>>radius;
			primitive *p=new sphere(pos,radius);
			pri.push_back(p);
		}
		else if(s=="OBJMESH"){
			char *file=new char[100];
			glm::vec3 pos;
			float scale;
			input>>s;
			input>>file;
			input>>s;
			input>>pos.x;input>>pos.y;input>>pos.z;
			input>>s;
			input>>scale;
			objmesh *obj=new objmesh(pos,scale);
			obj->read_from_file(file);
			primitive *p=obj;
			pri.push_back(p);
		}
		else if(s=="GPU_RENDER"){
			GPU_render=true;
		}
	}
	return true;
}

glm::vec3 fileInput::getCubeDim(){
	return cube_dim;
}

float fileInput::getSphereRadius(){
	return sphere_radius;
}

glm::vec3 fileInput::getParticleDim(){
	return particle_dim;
}

glm::vec3 fileInput::getParticleMinPos(){
	return particle_minPos;
}

glvuVec3f fileInput::getCameraLookAt(){
	return camera_lookat;
}

glvuVec3f fileInput::getCameraPos(){
	return camera_pos;
}

glvuVec3f fileInput::getCameraUp(){
	return camera_up;
}

float fileInput::getParticleDis(){
	return particle_dis;
}

int fileInput::getContainerType(){
	return container_type;
}

bool fileInput::isGPURender(){
	return GPU_render;
}