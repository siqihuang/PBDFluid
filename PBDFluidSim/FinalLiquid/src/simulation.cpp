#include "simulation.h"

simulation::simulation(){
	s=nullptr;
}

simulation::simulation(scene *s){
	this->s=s;
}

void simulation::updateOnCPU(){
	updateAcc();
	updateVel();
	updatePos();
	findNeighbor();
	for(int i=0;i<5;++i)
		PBDProjection();
	
	for(int i=0;i<s->particles.size();++i){
		s->particles[i].vel=(s->particles[i].newPos-s->particles[i].pos)/dt;
		s->particles[i].pos=s->particles[i].newPos;
	}
	
	s->cont->g->clearCellParticle();
}

void simulation::update(bool GPU_render){
	if(GPU_render) updateOnGPU();
	else updateOnCPU();
}

void simulation::updateAcc(){
	for(int i=0;i<s->particles.size();++i){
		s->particles[i].acc=gravityVector;
	}
}

void simulation::updateVel(){
	for(int i=0;i<s->particles.size();++i){
		s->particles[i].vel+=dt*s->particles[i].acc;
	}
}

void simulation::updatePos(){
	for(int i=0;i<s->particles.size();++i){
		s->particles[i].newPos=s->particles[i].pos+dt*s->particles[i].vel;
		s->cont->g->putParticle(s->particles[i],s->cont->minPos,i);
	}
	/*for(int i=0;i<20;++i){
		for(int j=0;j<20;++j){
			for(int k=0;k<20;++k){
				int index=i*400+j*20+k;
				if(s->cont->g->cells[i][j][k].cell_particles.size()!=0){
					//for(int m=0;m<s->cont->g->cells[i][j][k].cell_particles.size();++m)
						//cout<<s->cont->g->cells[i][j][k].cell_particles[m]<<","<<index<<endl;
					cout<<s->cont->g->cells[i][j][k].cell_particles.size()<<endl;
				}
			}
		}
	}*/
	//s->cont->g->print();
}

void simulation::collisionDetection(){
	for(int i=0;i<s->particles.size();++i){
		for(int j=0;j<s->pri.size();++j){
			s->pri[j]->staticIntersectionTest(s->particles[i]);
		}
		s->cont->staticIntersection(s->particles[i]);
	}
}

void simulation::findNeighbor(){
	int tmp=0;
	for(int i=0;i<s->particles.size();++i){
		s->particles[i].neighbor.clear();
		s->particles[i].neighbor=s->cont->g->neighbor_particles(i,s->particles[i].cellId);
		//int t=s->particles[i].neighbor.size();
		//tmp=max(tmp,t);
		//cout<<t<<","<<i<<endl;
	}
	//cout<<"@"<<endl;
	//cout<<tmp<<endl;
}

void simulation::PBDProjection(){
	for(int i=0;i<s->particles.size();++i){
		calculateLamda(s->particles[i]);
	}

	for(int i=0;i<s->particles.size();++i){
		glm::vec3 dp=glm::vec3(0,0,0);
		for(int j=0;j<s->particles[i].neighbor.size();++j){
			int index=s->particles[i].neighbor[j];
			glm::vec3 tmp=s->particles[i].newPos-s->particles[index].newPos;
			float lamdaCorrection=-1.0*0.1*pow(poly6Kernel(smooth_radius,s->particles[i].newPos-s->particles[index].newPos)/
				(poly6Kernel(smooth_radius,glm::vec3(0.1*smooth_radius,0,0))+Epsilon),4);
			lamdaCorrection=0;
			//if(abs(lamdaCorrection)>abs(s->particles[i].lamda+s->particles[index].lamda)*0.1) cout<<"inside"<<endl;
			dp+=(s->particles[i].lamda+s->particles[index].lamda+lamdaCorrection)*spikyKernel(smooth_radius,s->particles[i].newPos-s->particles[index].newPos)/rest_density;
		}
		s->particles[i].newPos=s->particles[i].newPos+dp;

		for(int j=0;j<s->pri.size();++j){
			s->pri[j]->staticIntersectionTest(s->particles[i]);
		}
		s->cont->staticIntersection(s->particles[i]);
		
	}
}

glm::vec3 simulation::spikyKernel(float h,glm::vec3 r){
	float R=glm::length(r);
	if(R>h) return glm::vec3(0,0,0);
	glm::vec3 result=(float)(15.0/(PI*pow(h,3))*pow(h-R,2))*r/(R+Epsilon);
	if(glm::length(result)>20) cout<<45.0/(PI*pow(h,3))*pow(h-R,2)<<endl;
	return (float)(15.0/(PI*pow(h,6))*pow(h-R,2))*r/(R+Epsilon);
}

float simulation::poly6Kernel(float h,glm::vec3 r){
	float R=glm::length(r);
	if(R>h) return 0;
	return 315.0/(64.0*PI*pow(h,9))*pow((h*h-glm::length(r)*glm::length(r)),2);
}

void simulation::calculateLamda(Particle &p){
	float SumGrediant=0.0;
	float ParticleGrediant=0.0;
	float density=0.0;
	for(int i=0;i<p.neighbor.size();++i){
		int index=p.neighbor[i];
		float tmp=-glm::length(spikyKernel(smooth_radius,p.newPos-s->particles[index].newPos)/rest_density);
		ParticleGrediant+=tmp;
		SumGrediant+=pow(tmp,2);
		density+=poly6Kernel(smooth_radius,p.newPos-s->particles[index].newPos);
	}
	SumGrediant += pow(ParticleGrediant,2);

	float densityConstriant=density/rest_density-1.0f;
	p.lamda=-1*densityConstriant/(SumGrediant+relaxation);
}