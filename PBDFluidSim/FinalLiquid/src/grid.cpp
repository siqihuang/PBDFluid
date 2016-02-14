#include "grid.h"


grid::grid(glm::vec3 boundingBox, float cell_size)
{
	grid_length = boundingBox[0];
	grid_height = boundingBox[1];
	grid_depth = boundingBox[2];
	this->cell_size = cell_size;
	gdim[0] = (int)ceil(grid_length/cell_size);
	gdim[1] = (int)ceil(grid_height/cell_size);
	gdim[2] = (int)ceil(grid_depth/cell_size);
	
	for(int i=0; i< gdim[0];i++){
		vector<vector<cell>> tmp1;
		cells.push_back(tmp1);
		for(int j=0;j<gdim[1];j++){
			vector<cell> tmp2;
			cells[i].push_back(tmp2);
			for(int k=0;k<gdim[2];k++){
				cell *tempC = new cell((i*grid_height*grid_depth+j*grid_depth+k),1);
				cells[i][j].push_back(*tempC);
			}
		}
	}
}


grid::~grid()
{
}

vector<int> grid::neighbor_particles(int particleID ,glm::vec3 cellID){
	vector<glm::vec3> temp_neighborcell; //vector stores all neighbor cells
	vector<int> ret; //return vectors which stores all neighbor particles of the input particle
	glm::vec3 temp_cell = cellID;
	for(int x=-1;x<2;x++){
		for(int y=-1;y<2;y++){
			for(int z=-1;z<2;z++){
				temp_cell[0]+=x;
				temp_cell[1]+=y;
				temp_cell[2]+=z;
				if(temp_cell[0]>=0&&temp_cell[0]<gdim[0]&&temp_cell[1]>=0&&temp_cell[1]<gdim[1]&&temp_cell[2]>=0&&temp_cell[2]<gdim[2])
					temp_neighborcell.push_back(temp_cell);
				temp_cell=cellID;
			}	
		}	
	}
	//cout<<temp_neighborcell.size()<<"!"<<endl;
	vector<int> curCell_particles;
	for(int i=0;i<temp_neighborcell.size();i++){
		glm::vec3 curNeighborCell = temp_neighborcell[i];
		//cout<<curNeighborCell<<endl;
		cell curCell = cells[curNeighborCell[0]][curNeighborCell[1]][curNeighborCell[2]];
		curCell_particles = curCell.get_particles();
		for(int j=0;j<curCell_particles.size();j++){
			if(particleID!=curCell_particles[j]){
				ret.push_back(curCell_particles[j]);
			}
		}
		curCell_particles.clear();
	}
	return ret;
}

float grid::get_length(){
	return grid_length;
}

float grid::get_height(){
	return grid_height;
}

float grid::get_depth(){
	return grid_depth;
}

void grid::print(){
	int tmp=0;
	for(int i=0;i<gdim[0];++i){
		for(int j=0;j<gdim[1];++j){
			for(int k=0;k<gdim[2];++k){
				//if(cells[i][j][k].cell_particles.size()>0)
				//cout<<cells[i][j][k].cell_particles.size()<<endl;
				int t=cells[i][j][k].cell_particles.size();
				tmp=max(tmp,t);
			}
		}
	}
	cout<<tmp<<endl;
}

void grid::putParticle(Particle &p,glm::vec3 offset,int particleId){
	glm::vec3 pos=p.newPos-offset;
	int x=(int)(pos[0]/cell_size);
	int y=(int)(pos[1]/cell_size);
	int z=(int)(pos[2]/cell_size);
	p.cellId=glm::vec3(x,y,z);
	//cout<<x<<","<<y<<","<<z<<endl;
	if(x<0||x>=gdim[0]||y<0||y>=gdim[1]||z<0||z>=gdim[2]) return;
	cells[x][y][z].add_particle(particleId);
}

void grid::clearCellParticle(){
	for(int i=0;i<gdim[0];++i){
		for(int j=0;j<gdim[1];++j){
			for(int k=0;k<gdim[2];++k){
				cells[i][j][k].clear_particles();
			}//k
		}//j
	}//i
}