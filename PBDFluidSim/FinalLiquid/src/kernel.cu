#include <thrust/execution_policy.h>
#include <thrust/random.h>
#include <thrust/device_vector.h>
#include <thrust/remove.h>
#include "kernel.h"

static int GparticleNum=0,GprimitiveNum=0;
static glm::vec3 GparticleDim=glm::vec3(0);
static GPUParticle *Gparticles,*dev_particles;
static GPUContainer *Gcontainer,*dev_container;
static GPUCell *Gcell,*dev_cell;
static GPUPrimitive *Gprimitive,*dev_primitive;

__host__ __device__ inline unsigned int utilhash(unsigned int a) {
    a = (a + 0x7ed55d16) + (a << 12);
    a = (a ^ 0xc761c23c) ^ (a >> 19);
    a = (a + 0x165667b1) + (a << 5);
    a = (a + 0xd3a2646c) ^ (a << 9);
    a = (a + 0xfd7046c5) + (a << 3);
    a = (a ^ 0xb55a4f09) ^ (a >> 16);
    return a;
}

__host__ __device__ thrust::default_random_engine random_engine(
        int iter, int index = 0, int depth = 0) {
    return thrust::default_random_engine(utilhash((index + 1) * iter) ^ utilhash(depth));
}

kdtree *initTree(kdtree *root){
	//postorder method to first get the left and right child on GPU Memory, then replace it with the memory on CPU, then copy the whole point to GPU
	if(root==nullptr) return nullptr;
	kdtree *dev_lc=initTree(root->lc);
	kdtree *dev_rc=initTree(root->rc);
	kdtree *tmp=new kdtree(root);
	tmp->lc=dev_lc;
	tmp->rc=dev_rc;
	kdtree *dev_root;
	cudaMalloc(&dev_root,sizeof(kdtree));
	cudaMemcpy(dev_root,tmp,sizeof(kdtree),cudaMemcpyHostToDevice);
	return dev_root;
}

__global__ void updateStateKernel(GPUParticle *particle,glm::vec3 acc,float dt,int N){
	int index=blockDim.x*blockIdx.x+threadIdx.x;
	if(index<N){
		particle[index].acc=acc;
		particle[index].vel+=dt*particle[index].acc;
		particle[index].newPos=particle[index].pos+dt*particle[index].vel;
		particle[index].index=0;
	}
}

__global__ void updateAccKernel(GPUParticle *particle,glm::vec3 acc,int N){
	int index=blockDim.x*blockIdx.x+threadIdx.x;
	if(index<N){
		particle[index].acc=acc;
	}
}

__global__ void updateVelKernel(GPUParticle *particle,float dt,int N){
	int index=blockDim.x*blockIdx.x+threadIdx.x;
	if(index<N){
		particle[index].vel+=dt*particle[index].acc;
	}
}

__global__ void updatePosKernel(GPUParticle *particle,float dt,int N){
	int index=blockDim.x*blockIdx.x+threadIdx.x;
	if(index<N){
		particle[index].newPos=particle[index].pos+dt*particle[index].vel;
	}
}

__global__ void clearCellParticleKernel(GPUContainer *container,int N){
	int index=blockDim.x*blockIdx.x+threadIdx.x;
	if(index<N){
		container->cell[index].index=0;
	}
}

__global__ void clearParticleNeighborKernel(GPUParticle *particle,int N){
	int index=blockDim.x*blockIdx.x+threadIdx.x;
	if(index<N){
		particle[index].index=0;
	}
}

__global__ void putParticleInCellKernel(GPUParticle *particle,GPUContainer *container,int cell_particle_limit,int N){
	int index=blockDim.x*blockIdx.x+threadIdx.x;
	if(index<N){
		glm::vec3 pos=particle[index].newPos-container->minPos;
		glm::vec3 cell_dim=container->cell_dim;
		pos.x=(int)(pos.x/container->cell_size);
		pos.y=(int)(pos.y/container->cell_size);
		pos.z=(int)(pos.z/container->cell_size);
		if(pos.x<0||pos.y<0||pos.z<0||pos.x>=cell_dim.x||pos.y>=cell_dim.y||pos.z>=cell_dim.z) return;
		int tmp=pos.x*cell_dim.y*cell_dim.z+pos.y*cell_dim.z+pos.z;
		//int value=container->cell[tmp].index;
		/*while(value!=atomicMax(&(container->cell[tmp].index),value)){
			value=container->cell[tmp].index;
		}*/
		/*
		this change is crucial
		*/
		int value=atomicAdd(&(container->cell[tmp].index),1);
		/*
		this change is crucial
		*/
		if(value>=cell_particle_limit){
			container->cell[tmp].index=cell_particle_limit;
			//return;
		}
		else{
			container->cell[tmp].particles[value]=index;
			particle[index].cellId=pos;
		}
	}
}

__global__ void findNeighborKernel(GPUParticle *particle,GPUContainer *container,int particle_neighbor_limit,int N){
	int index=blockDim.x*blockIdx.x+threadIdx.x;
	if(index<N){
		glm::vec3 cell_dim=container->cell_dim;
		int x=particle[index].cellId.x;
		int y=particle[index].cellId.y;
		int z=particle[index].cellId.z;
		for(int i=x-1;i<=x+1;++i){
			for(int j=y-1;j<=y+1;++j){
				for(int k=z-1;k<=z+1;++k){
					if(i<0||j<0||k<0||i>=cell_dim.x||j>=cell_dim.y||k>=cell_dim.z) continue;
					int tmp=cell_dim.y*cell_dim.z*i+cell_dim.z*j+k;
					for(int t=0;t<container->cell[tmp].index;++t){
						//if(container->cell[tmp].particles[t]!=index)
							particle[index].neighbor[particle[index].index++]=container->cell[tmp].particles[t];
						if(particle[index].index>=particle_neighbor_limit) return;
					}//t
				}//k
			}//j
		}//i
	}
}

__device__ glm::vec3 spikyKernelOnGPU(float h,float PI,glm::vec3 r){
	float R=glm::length(r);
	float Epsilon=1e-4;
	if(R>h) return glm::vec3(0,0,0);
	glm::vec3 result=(float)(15.0/(PI*pow(h,3))*pow(h-R,2))*r/(R+Epsilon);
	return (float)(15.0/(PI*pow(h,6))*pow(h-R,2))*r/(R+Epsilon);
}

__device__ float poly6KernelOnGPU(float h,float PI,glm::vec3 r){
	float R=glm::length(r);
	if(R>h) return 0;
	return 315.0/(64.0*PI*pow(h,9))*pow((h*h-glm::length(r)*glm::length(r)),2);
}

__device__ void cubeContainerIntersectionOnGPU(GPUParticle &p,float radius,glm::vec3 dimension,float dt,float rebound_rest,int index){
	glm::vec3 pos=p.newPos;
	glm::vec3 lower_bound=-dimension/2.0f,upper_bound=-lower_bound;
	for(int i=0;i<3;++i){
		//thrust::default_random_engine rng = random_engine(index, i, index+i);
        //thrust::uniform_real_distribution<float> u01(0, 1);
		if(pos[i]-radius-lower_bound[i]<=0){
			p.newPos[i]+=(radius+lower_bound[i]-pos[i]);
			//p.newPos[i]+=(0.04*1.0*u01(rng)-0.02);
			p.newPos[i]+=(0.04*1.0);
			p.vel[i]*=-rebound_rest;
			p.newPos[i]+=p.vel[i]*dt/2;
		}
		if(pos[i]+radius-upper_bound[i]>=0){
			p.newPos[i]-=(pos[i]+radius-upper_bound[i]);
			//p.newPos[i]+=(0.04*1.0*u01(rng)-0.02);
			p.newPos[i]-=(0.04*1.0);
			p.vel[i]*=-rebound_rest;
			p.newPos[i]+=p.vel[i]*dt/2;
		}
	}
}

__device__ void cubePrimitiveIntersectionOnGPU(GPUParticle &p,GPUPrimitive &primitive,float radius){
	glm::vec3 center = primitive.pos;
	glm::vec3 diff = p.newPos - center;
	glm::vec3 normal;
	float colision_epsilon=1e-3;
	float xcollide,ycollide,zcollide;
	glm::vec3 m_hf_dims=primitive.cubeDim/2.0f;
	xcollide=(fabs(diff.x)-m_hf_dims[0]-radius-colision_epsilon);
	ycollide=(fabs(diff.y)-m_hf_dims[1]-radius-colision_epsilon);
	zcollide=(fabs(diff.z)-m_hf_dims[2]-radius-colision_epsilon);
	if(xcollide<0&&ycollide<0&&zcollide<0){
		if(xcollide>=ycollide&&xcollide>=zcollide){
			if(diff.x>0){
				normal=glm::vec3(1,0,0);
			}
			else{
				normal=glm::vec3(-1,0,0);
			}
			p.newPos-=normal*xcollide;
		}
		else if(ycollide>=xcollide&&ycollide>=zcollide){
			if(diff.y>0){
				normal=glm::vec3(0,1,0);
			}
			else{
				normal=glm::vec3(0,-1,0);
			}
			p.newPos-=normal*ycollide;
		}
		else if(zcollide>=ycollide&&zcollide>=xcollide){
			if(diff.z>0){
				normal=glm::vec3(0,0,1);
			}
			else{
				normal=glm::vec3(0,0,-1);
			}
			p.newPos-=normal*zcollide;
		}
	}
}

__device__ void spherePrimitiveIntersectionOnGPU(GPUParticle &p,GPUPrimitive &primitive,float radius){
	glm::vec3 pos=p.newPos;
	float dis=glm::length(pos-primitive.pos);
	if(primitive.sphereRadius+radius>dis){
		dis=primitive.sphereRadius+radius-dis;
		glm::vec3 nor=glm::normalize(pos-primitive.pos);
		p.newPos+=nor*dis;

		//p.newPos+=(float)(0.04*1.0*rand()/(RAND_MAX+1)-0.02)*nor;
		//p.newPos-=glm::dot(p.vel,nor)*nor*dt/2.0f;

		glm::vec3 tmp=p.vel-glm::dot(p.vel,nor)*nor;
		p.vel=-glm::dot(p.vel,nor)*nor+tmp;
	}
}

__device__ bool insideBoxOnGPU(glm::vec3 pos,kdtree *tree){
	if(pos.x<=tree->xMax&&pos.x>=tree->xMin&&pos.y<=tree->yMax&&pos.y>=tree->yMin&&
		pos.z<=tree->zMax&&pos.z>=tree->zMin){
		return true;
	}
	else return false;
}

__device__ void getNearbyTrianglesOnGPU(glm::vec3 pos,kdtree *tree, int *list){
	int count=0,num=0,n=0;
	kdtree *kd[1000];
	kd[count++]=tree;
	while(count<1000&&n!=count&&num<180){
		kdtree *current=kd[n];
		if(insideBoxOnGPU(pos,current)){
			if(current->lc==nullptr&&current->rc==nullptr) list[num++]=current->index;
			else{
				kd[count++]=current->lc;
				if(count>=1000) break;
				kd[count++]=current->rc;
			}
		}
		n++;
	}
}

__device__ glm::vec3 getNormalOnGPU(glm::vec3 *m_positions,glm::vec3 *m_normals,int *m_indices, unsigned short TriangleIndex){
	glm::vec3 n1,n2,n3,v1,v2,v3,n,crossN,v12,v13;
	unsigned int index1,index2,index3;
	index1=m_indices[3*TriangleIndex];
	index2=m_indices[3*TriangleIndex+1];
	index3=m_indices[3*TriangleIndex+2];
	v1=m_positions[index1];v2=m_positions[index2];v3=m_positions[index3];
	n1=m_normals[index1];n2=m_normals[index2];n3=m_normals[index3];
	
	v12=v1-v2;v13=v1-v3;
	v12=glm::normalize(v12);v13=glm::normalize(v13);
	crossN=glm::cross(v12,v13);
	crossN=glm::normalize(crossN);
	
	n=(n1+n2+n3);
	n=glm::normalize(n);

	if(glm::dot(n,crossN)<0) return -crossN;
	else return crossN;
}

__device__ float getDistanceOnGPU(glm::vec3 *m_positions,glm::vec3 *m_normals,int *m_indices,glm::vec3 p,unsigned short TriangleIndex){
	float dis,k,x;
	unsigned short index;
	index=m_indices[3*TriangleIndex];
	
	glm::vec3 normal=getNormalOnGPU(m_positions,m_normals,m_indices,TriangleIndex);
	
	glm::vec3 d=p-m_positions[index];
	x=(normal.x*d.x+normal.y*d.y+normal.z*d.z);
	return x;
}

__device__ void objmeshPrimitiveIntersectionOnGPU(GPUParticle &p,GPUPrimitive &primitive,float radius){
	float minDis=-1e7;
	float COLLISION_EPSILON=1e-2+radius;
	bool inCollision=false;
	glm::vec3 pos=p.newPos;
	glm::vec3 normal(0);
	int list[180];
	for(int i=0;i<180;++i) list[i]=-1;
	getNearbyTrianglesOnGPU(pos,primitive.mesh,list);
	pos-=primitive.pos;

	for(int i=0;i<180;i++){
		if(list[i]==-1) break;
		float tmp=getDistanceOnGPU(primitive.objVertex,primitive.objNormal,primitive.objIndices,pos,list[i])-COLLISION_EPSILON;
		if(tmp<0&&tmp>minDis&&tmp>-0.5){
			glm::vec3 n=getNormalOnGPU(primitive.objVertex,primitive.objNormal,primitive.objIndices,list[i]);
			normal=n;
			minDis=tmp;
			inCollision=true;
		}
	}
	if(inCollision){
		p.newPos-=normal*minDis;
		p.pos=p.newPos;
		glm::vec3 tmp=p.vel-glm::dot(p.vel,normal)*normal;
		p.vel=(tmp-glm::dot(p.vel,normal)*normal)*0.6f;
		p.newPos+=0.01f*p.vel/4.0f;
	}
	
}

__global__ void primitiveCollisionKernel(GPUParticle *particle,GPUPrimitive *primitive,int primitiveNum,float radius,int N){
	int index=blockDim.x*blockIdx.x+threadIdx.x;
	if(index<N){
		for(int i=0;i<primitiveNum;++i){
			if(primitive[i].primitive_type==0){//cube
				cubePrimitiveIntersectionOnGPU(particle[index],primitive[i],radius);
			}
			else if(primitive[i].primitive_type==1){//sphere
				spherePrimitiveIntersectionOnGPU(particle[index],primitive[i],radius);
			}
			else if(primitive[i].primitive_type==2){//objmesh
				objmeshPrimitiveIntersectionOnGPU(particle[index],primitive[i],radius);
			}
		}
	}
}

__global__ void calculateLamdaKernel(GPUParticle *particle,float rest_density,float smooth_radius,float relaxation,float PI,int N){
	int index=blockDim.x*blockIdx.x+threadIdx.x;
	if(index<N){
		float SumGrediant=0.0;
		float ParticleGrediant=0.0;
		float density=0.0;
		for(int i=0;i<particle[index].index;++i){
			int index1=particle[index].neighbor[i];
			float tmp=-glm::length(spikyKernelOnGPU(smooth_radius,PI,particle[index].newPos-particle[index1].newPos)/rest_density);
			ParticleGrediant+=tmp;
			SumGrediant+=pow(tmp,2);
			density+=poly6KernelOnGPU(smooth_radius,PI,particle[index].newPos-particle[index1].newPos);
		}
		SumGrediant += pow(ParticleGrediant,2);
	
		float densityConstriant=density/rest_density-1.0f;
		particle[index].lamda=-1*densityConstriant/(SumGrediant+relaxation);
	}
}

__global__ void PBDProjectionKernel(GPUParticle *particle,float rest_density,float smooth_radius,float PI,int N){
	int index=blockDim.x*blockIdx.x+threadIdx.x;
	if(index<N){
		glm::vec3 dp=glm::vec3(0,0,0);
		for(int j=0;j<particle[index].index;++j){
			int index1=particle[index].neighbor[j];
			glm::vec3 tmp=particle[index].newPos-particle[index1].newPos;
			dp+=(particle[index].lamda+particle[index1].lamda)*spikyKernelOnGPU(smooth_radius,PI,particle[index].newPos-particle[index1].newPos)/rest_density;
		}
		particle[index].newPos=particle[index].newPos+dp;
			
		/*for(int j=0;j<s->pri.size();++j){
			s->pri[j]->staticIntersectionTest(s->particles[i]);
		}
		s->cont->staticIntersection(s->particles[i]);*/
	}
}

__global__ void containerCollisionKernel(GPUParticle *particle,GPUContainer *container,float particle_radius,float dt,float rebound_rest,int N){
	int index=blockDim.x*blockIdx.x+threadIdx.x;
	if(index<N){
		if(container->container_type==0){//CUBE
			cubeContainerIntersectionOnGPU(particle[index],particle_radius,container->cubeDim,dt,rebound_rest,index);
		}
	}
}

__global__ void calculateNewVelocityKernel(GPUParticle *particle,float dt,int N){
	int index=blockDim.x*blockIdx.x+threadIdx.x;
	if(index<N){
		particle[index].vel=(particle[index].newPos-particle[index].pos)/dt;
		particle[index].pos=particle[index].newPos;
	}
}

void initParticleOnGPU(glm::vec3 particleDim,glm::vec3 particleMinPos,float particleDis){
	GparticleDim=particleDim;
	GparticleNum=GparticleDim[0]*GparticleDim[1]*GparticleDim[2];
	Gparticles=new GPUParticle[GparticleNum];
	for(int i=0;i<GparticleDim.x;++i){
		for(int j=0;j<GparticleDim.y;++j){
			for(int k=0;k<GparticleDim.z;++k){
				int index=i*GparticleDim.y*GparticleDim.z+j*GparticleDim.z+k;
				Gparticles[index].acc=glm::vec3(0);
				Gparticles[index].vel=glm::vec3(0);
				Gparticles[index].pos=particleMinPos+particleDis*glm::vec3(i,j,k);
				Gparticles[index].newPos=glm::vec3(0);
				Gparticles[index].cellId=glm::vec3(0);
				Gparticles[index].lamda=0;
				Gparticles[index].mass=0;
				Gparticles[index].index=0;
			}//k
		}//j
	}
	cudaMalloc(&dev_particles,GparticleNum*sizeof(GPUParticle));
	cudaMemcpy(dev_particles,Gparticles,GparticleNum*sizeof(GPUParticle),cudaMemcpyHostToDevice);
}

void initCubeContainerOnGPU(glm::vec3 dimension,float cell_size){
	Gcontainer=new GPUContainer();
	Gcontainer->container_type=CONTAINER::CONTAINER_CUBE;
	Gcontainer->cell_size=cell_size;
	Gcontainer->cubeDim=dimension;
	Gcontainer->boundingBox=dimension;
	Gcontainer->minPos=-dimension/2.0f;
	int x=(int)ceil(Gcontainer->boundingBox.x/cell_size);
	int y=(int)ceil(Gcontainer->boundingBox.y/cell_size);
	int z=(int)ceil(Gcontainer->boundingBox.z/cell_size);
	Gcontainer->cell_dim=glm::vec3(x,y,z);

	initCellOnGPU(Gcontainer->boundingBox,Gcontainer->cell_dim,cell_size);
	Gcontainer->cell=dev_cell;

	cudaMalloc(&dev_container,sizeof(GPUContainer));
	cudaMemcpy(dev_container,Gcontainer,sizeof(GPUContainer),cudaMemcpyHostToDevice);
}

void initSphereContainerOnGPU(float radius,float cell_size){
	Gcontainer=new GPUContainer();
	Gcontainer->container_type=CONTAINER::CONTAINER_SPHERE;
	Gcontainer->cell_size=cell_size;
	Gcontainer->sphereRadius=radius;
	Gcontainer->boundingBox=2.0f*glm::vec3(radius);
	Gcontainer->minPos=glm::vec3(-radius);
	int x=(int)ceil(Gcontainer->boundingBox.x/cell_size);
	int y=(int)ceil(Gcontainer->boundingBox.y/cell_size);
	int z=(int)ceil(Gcontainer->boundingBox.z/cell_size);
	Gcontainer->cell_dim=glm::vec3(x,y,z);

	initCellOnGPU(Gcontainer->boundingBox,Gcontainer->cell_dim,cell_size);
	Gcontainer->cell=dev_cell;

	cudaMalloc(&dev_container,sizeof(GPUContainer));
	cudaMemcpy(dev_container,Gcontainer,sizeof(GPUContainer),cudaMemcpyHostToDevice);
}

void initCellOnGPU(glm::vec3 boundingBox,glm::vec3 cell_dim,float cell_size){
	int x=cell_dim.x;
	int y=cell_dim.y;
	int z=cell_dim.z;
	Gcell=new GPUCell[x*y*z];
	for(int i=0;i<x;++i){
		for(int j=0;j<y;++j){
			for(int k=0;k<z;++k){
				int index=i*y*z+j*z+k;
				Gcell[index].pos=-boundingBox/2.0f+glm::vec3(i,j,k)*cell_size;
				Gcell[index].index=0;
			}//k
		}//j
	}//i
	cudaMalloc(&dev_cell,x*y*z*sizeof(GPUCell));
	cudaMemcpy(dev_cell,Gcell,x*y*z*sizeof(GPUCell),cudaMemcpyHostToDevice);
}

void initPrimitiveOnGPU(GPUPrimitive *primitive,int num){
	Gprimitive=primitive;
	GprimitiveNum=num;
	cudaMalloc(&dev_primitive,num*sizeof(GPUPrimitive));
	cudaMemcpy(dev_primitive,Gprimitive,num*sizeof(GPUPrimitive),cudaMemcpyHostToDevice);
}

void updateOnGPU(){
	//updateAccOnGPU();
	//updateVelOnGPU();
	//updatePosOnGPU();
	//clearParticleNeighborOnGPU();
	updateStateOnGPU();

	clearCellParticleOnGPU();
	putParticleInCellOnGPU();
	findNeighborOnGPU();
	for(int i=0;i<4;++i){
		calculateLamdaOnGPU();
		PBDProjectionOnGPU();
		primitiveCollisionOnGPU();
		containerCollisionOnGPU();
		//cudaMemcpy(Gparticles,dev_particles,GparticleNum*sizeof(GPUParticle),cudaMemcpyDeviceToHost);
		//cout<<"newPos "<<Gparticles[0].newPos.x<<","<<Gparticles[0].newPos.y<<","<<Gparticles[0].newPos.z<<endl;
		//cout<<"vel "<<Gparticles[0].vel.x<<","<<Gparticles[0].vel.y<<","<<Gparticles[0].vel.z<<endl;
	}
	calculateNewVelOnGPU();
	cudaMemcpy(Gparticles,dev_particles,GparticleNum*sizeof(GPUParticle),cudaMemcpyDeviceToHost);
}

void updateAccOnGPU(){
	updateAccKernel<<<(GparticleNum+255)/256,256>>>(dev_particles,gravityVector,GparticleNum);
}

void updateVelOnGPU(){
	updateVelKernel<<<(GparticleNum+255)/256,256>>>(dev_particles,dt,GparticleNum);
}

void updatePosOnGPU(){
	updatePosKernel<<<(GparticleNum+255)/256,256>>>(dev_particles,dt,GparticleNum);
}

void putParticleInCellOnGPU(){
	putParticleInCellKernel<<<(GparticleNum+255)/256,256>>>(dev_particles,dev_container,GPU_cell_particle_limit,GparticleNum);
	//cudaMemcpy(Gcontainer,dev_container,sizeof(GPUContainer),cudaMemcpyDeviceToHost);
	//cudaMemcpy(Gcell,dev_cell,20*20*20*sizeof(GPUCell),cudaMemcpyDeviceToHost);
	/*for(int i=0;i<20*20*20;++i){
		if(Gcell[i].index!=0){
			//for(int j=0;j<Gcell[i].index;++j)
				//cout<<Gcell[i].particles[j]<<","<<i<<endl;
			cout<<Gcell[i].index<<endl;
		}
	}*/
}

void clearCellParticleOnGPU(){
	int num=Gcontainer->cell_dim.x*Gcontainer->cell_dim.y*Gcontainer->cell_dim.z;
	clearCellParticleKernel<<<(num+255)/256,256>>>(dev_container,num);
}

void clearParticleNeighborOnGPU(){
	clearParticleNeighborKernel<<<(GparticleNum+255)/256,256>>>(dev_particles,GparticleNum);
}

void findNeighborOnGPU(){
	findNeighborKernel<<<(GparticleNum+255)/256,256>>>(dev_particles,dev_container,GPU_particle_neighbor_limit,GparticleNum);
	//cudaMemcpy(Gparticles,dev_particles,GparticleNum*sizeof(GPUParticle),cudaMemcpyDeviceToHost);
	/*for(int i=0;i<GparticleNum;++i){
		cout<<Gparticles[i].index<<","<<i<<endl;
	}
	cout<<"@"<<endl;*/
}

void calculateLamdaOnGPU(){
	calculateLamdaKernel<<<(GparticleNum+255)/256,256>>>(dev_particles,rest_density,smooth_radius,relaxation,PI,GparticleNum);
}

void PBDProjectionOnGPU(){
	PBDProjectionKernel<<<(GparticleNum+255)/256,256>>>(dev_particles,rest_density,smooth_radius,PI,GparticleNum);
}

void calculateNewVelOnGPU(){
	calculateNewVelocityKernel<<<(GparticleNum+255)/256,256>>>(dev_particles,dt,GparticleNum);
}

void containerCollisionOnGPU(){
	containerCollisionKernel<<<(GparticleNum+255)/256,256>>>(dev_particles,dev_container,particle_radius,dt,rebound_rest,GparticleNum);
}

void primitiveCollisionOnGPU(){
	/*glm::vec3 *tmp,*dev_tmp;
	float *tmp1,*dev_tmp1;
	tmp=new glm::vec3[1];
	tmp1=new float[1];
	cudaMalloc(&dev_tmp,sizeof(glm::vec3));
	cudaMalloc(&dev_tmp1,sizeof(float));*/
	primitiveCollisionKernel<<<(GparticleNum+255)/256,256>>>(dev_particles,dev_primitive,GprimitiveNum,particle_radius,GparticleNum);
	/*cudaMemcpy(tmp,dev_tmp,sizeof(glm::vec3),cudaMemcpyDeviceToHost);
	cudaMemcpy(tmp1,dev_tmp1,sizeof(float),cudaMemcpyDeviceToHost);
	cout<<tmp[0].x<<","<<tmp[0].y<<","<<tmp[0].z<<endl;
	cout<<tmp1[0]<<endl;
	delete(tmp);
	delete(tmp1);*/
}

glm::vec3 getPos(int n){
	return Gparticles[n].pos;
}

void updateStateOnGPU(){
	updateStateKernel<<<(GparticleNum+255)/256,256>>>(dev_particles,gravityVector,dt,GparticleNum);
}