#include "primitive.h"

primitive::primitive(){
}

primitive::primitive(glm::vec3 pos){
	this->pos=pos;
}

void primitive::init(){
}

void primitive::staticIntersectionTest(Particle &p){
}

/*
borrowed from tiantian Liu's deformable_body_sim
*/
void primitive::draw(const VBO& vbos){	
	glBindBuffer(GL_ARRAY_BUFFER, vbos.m_vbo);
    glBufferData(GL_ARRAY_BUFFER, 3 * m_positions.size() * sizeof(float), &m_positions[0], GL_STREAM_DRAW);
    // color
    glBindBuffer(GL_ARRAY_BUFFER, vbos.m_cbo);
    glBufferData(GL_ARRAY_BUFFER, 3 * m_colors.size() * sizeof(float), &m_colors[0], GL_STREAM_DRAW);

    // normal
    glBindBuffer(GL_ARRAY_BUFFER, vbos.m_nbo);
    glBufferData(GL_ARRAY_BUFFER, 3 * m_normals.size() * sizeof(float), &m_normals[0], GL_STREAM_DRAW);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbos.m_ibo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, m_indices.size() * sizeof(unsigned short), &m_indices[0], GL_STATIC_DRAW);

    glEnableVertexAttribArray(0);
    glEnableVertexAttribArray(1);
    glEnableVertexAttribArray(2);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbos.m_ibo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, m_indices.size() * sizeof(unsigned short), &m_indices[0], GL_STATIC_DRAW);

    glEnableVertexAttribArray(0);
    glEnableVertexAttribArray(1);
    glEnableVertexAttribArray(2);

    glBindBuffer(GL_ARRAY_BUFFER, vbos.m_vbo);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

    glBindBuffer(GL_ARRAY_BUFFER, vbos.m_cbo);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, 0);

    glBindBuffer(GL_ARRAY_BUFFER, vbos.m_nbo);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 0, 0);

    glm::mat4 transformation(1.0f);
    transformation = glm::translate(transformation, pos);

    glUniformMatrix4fv(vbos.m_uniform_transformation, 1, false, &transformation[0][0]);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbos.m_ibo);
    glDrawElements(GL_TRIANGLES, m_indices.size(), GL_UNSIGNED_SHORT, 0);//GL_UNSIGNED_INT

    glDisableVertexAttribArray(0);
    glDisableVertexAttribArray(1);
    glDisableVertexAttribArray(2);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	//cout<<pos.x<<","<<pos.y<<","<<pos.z<<endl;
}

GPUPrimitive primitive::getGPUPrimitive(){
	GPUPrimitive p;
	return p;
}

cube::cube(){}

cube::cube(glm::vec3 pos,glm::vec3 dimension):primitive(pos){
	this->dimension=dimension;
	this->m_hf_dims=dimension/2.0f;
	this->lower_bound=pos-m_hf_dims;
	this->upper_bound=pos+m_hf_dims;
	this->primitive_type=PRIMITIVE::CUBE;
	this->init();
}

void cube::init(){
	m_positions.clear();
    m_colors.clear();
    m_normals.clear();
    m_indices.clear();

    glm::vec3 mat_color(0.6);

    // front face 012, 321
    m_positions.push_back(glm::vec3(-m_hf_dims.x, -m_hf_dims.y, +m_hf_dims.z));
    m_colors.push_back(mat_color);
    m_normals.push_back(glm::vec3(0,0,1));
    m_positions.push_back(glm::vec3(+m_hf_dims.x, -m_hf_dims.y, +m_hf_dims.z));
    m_colors.push_back(mat_color);
    m_normals.push_back(glm::vec3(0,0,1));
    m_positions.push_back(glm::vec3(-m_hf_dims.x, +m_hf_dims.y, +m_hf_dims.z));
    m_colors.push_back(mat_color);
    m_normals.push_back(glm::vec3(0,0,1));
    m_positions.push_back(glm::vec3(+m_hf_dims.x, +m_hf_dims.y, +m_hf_dims.z));
    m_colors.push_back(mat_color);
    m_normals.push_back(glm::vec3(0,0,1));
    m_indices.push_back(0);
    m_indices.push_back(1);
    m_indices.push_back(2);
    m_indices.push_back(3);
    m_indices.push_back(2);
    m_indices.push_back(1);

    // back face 654 567
    m_positions.push_back(glm::vec3(-m_hf_dims.x, -m_hf_dims.y, -m_hf_dims.z));
    m_colors.push_back(mat_color);
    m_normals.push_back(glm::vec3(0,0,-1));
    m_positions.push_back(glm::vec3(+m_hf_dims.x, -m_hf_dims.y, -m_hf_dims.z));
    m_colors.push_back(mat_color);
    m_normals.push_back(glm::vec3(0,0,-1));
    m_positions.push_back(glm::vec3(-m_hf_dims.x, +m_hf_dims.y, -m_hf_dims.z));
    m_colors.push_back(mat_color);
    m_normals.push_back(glm::vec3(0,0,-1));
    m_positions.push_back(glm::vec3(+m_hf_dims.x, +m_hf_dims.y, -m_hf_dims.z));
    m_colors.push_back(mat_color);
    m_normals.push_back(glm::vec3(0,0,-1));
    m_indices.push_back(6);
    m_indices.push_back(5);
    m_indices.push_back(4);
    m_indices.push_back(5);
    m_indices.push_back(6);
    m_indices.push_back(7);

    // right face 8 9 10, 11 10 9
    m_positions.push_back(glm::vec3(+m_hf_dims.x, -m_hf_dims.y, +m_hf_dims.z));
    m_colors.push_back(mat_color);
    m_normals.push_back(glm::vec3(1,0,0));
    m_positions.push_back(glm::vec3(+m_hf_dims.x, -m_hf_dims.y, -m_hf_dims.z));
    m_colors.push_back(mat_color);
    m_normals.push_back(glm::vec3(1,0,0));
    m_positions.push_back(glm::vec3(+m_hf_dims.x, +m_hf_dims.y, +m_hf_dims.z));
    m_colors.push_back(mat_color);
    m_normals.push_back(glm::vec3(1,0,0));
    m_positions.push_back(glm::vec3(+m_hf_dims.x, +m_hf_dims.y, -m_hf_dims.z));
    m_colors.push_back(mat_color);
    m_normals.push_back(glm::vec3(1,0,0));
    m_indices.push_back(8);
    m_indices.push_back(9);
    m_indices.push_back(10);
    m_indices.push_back(11);
    m_indices.push_back(10);
    m_indices.push_back(9);

    // left face 14 13 12, 13 14 15
    m_positions.push_back(glm::vec3(-m_hf_dims.x, -m_hf_dims.y, +m_hf_dims.z));
    m_colors.push_back(mat_color);
    m_normals.push_back(glm::vec3(-1,0,0));
    m_positions.push_back(glm::vec3(-m_hf_dims.x, -m_hf_dims.y, -m_hf_dims.z));
    m_colors.push_back(mat_color);
    m_normals.push_back(glm::vec3(-1,0,0));
    m_positions.push_back(glm::vec3(-m_hf_dims.x, +m_hf_dims.y, +m_hf_dims.z));
    m_colors.push_back(mat_color);
    m_normals.push_back(glm::vec3(-1,0,0));
    m_positions.push_back(glm::vec3(-m_hf_dims.x, +m_hf_dims.y, -m_hf_dims.z));
    m_colors.push_back(mat_color);
    m_normals.push_back(glm::vec3(-1,0,0));
    m_indices.push_back(14);
    m_indices.push_back(13);
    m_indices.push_back(12);
    m_indices.push_back(13);
    m_indices.push_back(14);
    m_indices.push_back(15);

    // top face 16 17 18, 19 18 17
    m_positions.push_back(glm::vec3(-m_hf_dims.x, +m_hf_dims.y, +m_hf_dims.z));
    m_colors.push_back(mat_color);
    m_normals.push_back(glm::vec3(0,1,0));
    m_positions.push_back(glm::vec3(+m_hf_dims.x, +m_hf_dims.y, +m_hf_dims.z));
    m_colors.push_back(mat_color);
    m_normals.push_back(glm::vec3(0,1,0));
    m_positions.push_back(glm::vec3(-m_hf_dims.x, +m_hf_dims.y, -m_hf_dims.z));
    m_colors.push_back(mat_color);
    m_normals.push_back(glm::vec3(0,1,0));
    m_positions.push_back(glm::vec3(+m_hf_dims.x, +m_hf_dims.y, -m_hf_dims.z));
    m_colors.push_back(mat_color);
    m_normals.push_back(glm::vec3(0,1,0));
    m_indices.push_back(16);
    m_indices.push_back(17);
    m_indices.push_back(18);
    m_indices.push_back(19);
    m_indices.push_back(18);
    m_indices.push_back(17);

    // bottom face 22 21 20, 21 22 23
    m_positions.push_back(glm::vec3(-m_hf_dims.x, -m_hf_dims.y, +m_hf_dims.z));
    m_colors.push_back(mat_color);
    m_normals.push_back(glm::vec3(0,-1,0));
    m_positions.push_back(glm::vec3(+m_hf_dims.x, -m_hf_dims.y, +m_hf_dims.z));
    m_colors.push_back(mat_color);
    m_normals.push_back(glm::vec3(0,-1,0));
    m_positions.push_back(glm::vec3(-m_hf_dims.x, -m_hf_dims.y, -m_hf_dims.z));
    m_colors.push_back(mat_color);
    m_normals.push_back(glm::vec3(0,-1,0));
    m_positions.push_back(glm::vec3(+m_hf_dims.x, -m_hf_dims.y, -m_hf_dims.z));
    m_colors.push_back(mat_color);
    m_normals.push_back(glm::vec3(0,-1,0));
    m_indices.push_back(22);
    m_indices.push_back(21);
    m_indices.push_back(20);
    m_indices.push_back(21);
    m_indices.push_back(22);
    m_indices.push_back(23);
}

void cube::staticIntersectionTest(Particle &p){
	glm::vec3 center = pos;
	glm::vec3 diff = p.newPos - center;
	glm::vec3 normal;
	float radius=p.getRadius();
	float xcollide,ycollide,zcollide;
	xcollide=(fabs(diff.x)-this->m_hf_dims[0]-radius-colision_epsilon);
	ycollide=(fabs(diff.y)-this->m_hf_dims[1]-radius-colision_epsilon);
	zcollide=(fabs(diff.z)-this->m_hf_dims[2]-radius-colision_epsilon);
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

GPUPrimitive cube::getGPUPrimitive(){
	GPUPrimitive p;
	p.primitive_type=PRIMITIVE::CUBE;
	p.cubeDim=dimension;
	p.pos=pos;
	return p;
}

sphere::sphere(){}

sphere::sphere(glm::vec3 pos,float radius):primitive(pos){
	this->radius=radius;
	this->primitive_type=PRIMITIVE::SPHERE;
	this->init();
}

void sphere::init(){
    m_positions.clear();
    m_colors.clear();
    m_normals.clear();
    m_indices.clear();

    glm::vec3 mat_color(0.6);
    unsigned int slice = 24, stack = 10;

    glm::vec3 tnormal(0.0, 1.0, 0.0), tpos;
    tpos = radius * tnormal;

    m_positions.push_back(tpos);
    m_normals.push_back(tnormal);
    m_colors.push_back(mat_color);

    float theta_z, theta_y, sin_z;
    float delta_y = 360.0 / slice, delta_z = 180.0 / stack;
    //loop over the sphere
    for(theta_z = delta_z; theta_z < 179.99; theta_z += delta_z)
    {
        for(theta_y = 0.0; theta_y < 359.99; theta_y += delta_y)
        {
            sin_z = sin(glm::radians(theta_z));
            
            tnormal.x = sin_z * cos(glm::radians(theta_y));
            tnormal.y = cos(glm::radians(theta_z));
            tnormal.z = -sin_z * sin(glm::radians(theta_y));

            tpos = radius * tnormal;

            m_positions.push_back(tpos);
            m_normals.push_back(tnormal);
            m_colors.push_back(mat_color);
        }
    }
    tnormal = glm::vec3(0.0, -1.0, 0.0);
    tpos = radius * tnormal;

    m_positions.push_back(tpos);
    m_normals.push_back(tnormal);
    m_colors.push_back(mat_color);

    //indices
    unsigned int j = 0, k = 0;
    for(j = 0; j < slice - 1; ++j)
    {
        m_indices.push_back(0);
        m_indices.push_back(j + 1);
        m_indices.push_back(j + 2);
    }
    m_indices.push_back(0);
    m_indices.push_back(slice);
    m_indices.push_back(1);

    for(j = 0; j < stack - 2; ++j)
    {
        for(k = 1 + slice * j; k < slice * (j + 1); ++k)
        {
            m_indices.push_back(k);
            m_indices.push_back(k + slice);
            m_indices.push_back(k + slice + 1);

            m_indices.push_back(k);
            m_indices.push_back(k + slice + 1);
            m_indices.push_back(k + 1);
        }
        m_indices.push_back(k);
        m_indices.push_back(k + slice);
        m_indices.push_back(k + 1);

        m_indices.push_back(k);
        m_indices.push_back(k + 1);
        m_indices.push_back(k + 1 - slice);
    }

    unsigned int bottom_id = (stack - 1) * slice + 1;
    unsigned int offset = bottom_id - slice;
    for(j = 0; j < slice - 1; ++j)
    {
        m_indices.push_back(j + offset);
        m_indices.push_back(bottom_id);
        m_indices.push_back(j + offset + 1);
    }
    m_indices.push_back(bottom_id - 1);
    m_indices.push_back(bottom_id);
    m_indices.push_back(offset);

    if(m_indices.size() != 6 * (stack - 1) * slice)
        printf("indices number not correct!\n");
}

void sphere::staticIntersectionTest(Particle &p){
	glm::vec3 pos=p.newPos;
	float radius=p.getRadius();
	float dis=glm::length(pos-this->pos);
	if(this->radius+radius>dis){
		dis=this->radius+radius-dis;
		glm::vec3 nor=glm::normalize(pos-this->pos);
		p.newPos+=nor*dis;

		//p.newPos+=(float)(0.04*1.0*rand()/(RAND_MAX+1)-0.02)*nor;
		p.newPos-=glm::dot(p.vel,nor)*nor*dt/2.0f;

		glm::vec3 tmp=p.vel-glm::dot(p.vel,nor)*nor;
		p.vel=-glm::dot(p.vel,nor)*nor*rebound_rest+rebound_rest*tmp;
	}
}

GPUPrimitive sphere::getGPUPrimitive(){
	GPUPrimitive p;
	p.primitive_type=PRIMITIVE::SPHERE;
	p.sphereRadius=radius;
	p.pos=pos;
	return p;
}

objmesh::objmesh(){}

objmesh::objmesh(glm::vec3 pos,float m_scaling):primitive(pos){
	this->m_scaling=m_scaling;
	this->primitive_type=PRIMITIVE::OBJMESH;
	this->init();
}

void objmesh::init(){

}

void objmesh::staticIntersectionTest(Particle &p){
}

GPUPrimitive objmesh::getGPUPrimitive(){
	GPUPrimitive p;
	p.primitive_type=PRIMITIVE::OBJMESH;
	glm::vec3 *dev_objVertex,*dev_objNormal,*objVertex,*objNormal;
	int *dev_objIndices,*objIndices;
	objVertex=new glm::vec3[m_positions.size()];
	objNormal=new glm::vec3[m_normals.size()];
	objIndices=new int[m_indices.size()];
	cudaMalloc(&dev_objVertex,m_positions.size()*sizeof(glm::vec3));
	cudaMalloc(&dev_objNormal,m_normals.size()*sizeof(glm::vec3));
	cudaMalloc(&dev_objIndices,m_indices.size()*sizeof(int));
	for(int i=0;i<m_positions.size();++i){
		objVertex[i]=m_positions[i];
	}
	for(int i=0;i<m_normals.size();++i){
		objNormal[i]=m_normals[i];
	}
	for(int i=0;i<m_indices.size();++i){
		objIndices[i]=m_indices[i];
	}
	cudaMemcpy(dev_objVertex,objVertex,m_positions.size()*sizeof(glm::vec3),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_objNormal,objNormal,m_normals.size()*sizeof(glm::vec3),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_objIndices,objIndices,m_indices.size()*sizeof(int),cudaMemcpyHostToDevice);
	p.objVertex=dev_objVertex;
	p.objNormal=dev_objNormal;
	p.objIndices=dev_objIndices;

	p.mesh=initTree(&tree);
	p.pos=pos;
	return p;
}

void objmesh::read_from_file(char* filename)
{
    m_positions.clear();
    m_colors.clear();
    m_normals.clear();
    m_indices.clear();
    glm::vec3 mat_color(0.6); 
    // vertices and color
	
    std::ifstream infile(filename);
    if(!infile.good())
    {
        printf("Error in loading file %s\n", filename);
        exit(0);
    }
    char buffer[256];
    unsigned int ip0, ip1, ip2;
    unsigned int n0, n1, n2;
    glm::vec3 pos;

    while(!infile.getline(buffer,255).eof())
    {
        buffer[255] = '\0';
        if(buffer[0] == 'v' && (buffer[1] == ' ' || buffer[1] == 32))
        {
            if(sscanf_s(buffer, "v %f %f %f", &pos.x, &pos.y, &pos.z) == 3)
            {
                pos = m_scaling * pos;
                m_positions.push_back(pos);
            }
            else
            {
                printf("Vertex is not in desired format.\n");
                exit(0);
            }
        }
        else if (buffer[0] == 'v' && buffer[1] == 'n' && (buffer[2] == ' ' || buffer[2] == 32))
        {
            // load normals from obj file.
        }
        else if (buffer[0] == 'f' && (buffer[1] == ' ' || buffer[1] == 32))
        {
            if(sscanf_s(buffer, "f %u %u %u", &ip0, &ip1, &ip2) == 3)
            {
                m_indices.push_back(--ip0);
                m_indices.push_back(--ip1);
                m_indices.push_back(--ip2);
            }
            else if(sscanf_s(buffer, "f %u//%u %u//%u %u//%u", &ip0, &n0, &ip1, &n1, &ip2, &n2) == 6)
            {
                m_indices.push_back(--ip0);
                m_indices.push_back(--ip1);
                m_indices.push_back(--ip2);
            }
            else if(sscanf_s(buffer, "f %u/%u %u/%u %u/%u", &ip0, &n0, &ip1, &n1, &ip2, &n2) == 6)
            {
                m_indices.push_back(--ip0);
                m_indices.push_back(--ip1);
                m_indices.push_back(--ip2);
            }
            else
            {
                printf("Triangle indices is not in desired format.\n");
                exit(0);
            }
        }
    }
    // normals

    unsigned int id, size;
    bool vert_norm = (m_normals.size() != m_positions.size());
    if(vert_norm)
        m_normals.resize(m_positions.size(), glm::vec3(0.0f));

    size = m_indices.size();
    glm::uvec3 triangle;
    glm::vec3 p0, p1, p2;
    glm::vec3 norm;
    float phi0, phi1, phi2;
    float pi = glm::radians(180.0f);
    for(id = 0; id < size; id+=3)
    {
        triangle = glm::uvec3(m_indices[id], m_indices[id+1], m_indices[id+2]);
        p0 = m_positions[triangle.x];
        p1 = m_positions[triangle.y];
        p2 = m_positions[triangle.z];
        norm = glm::normalize(glm::cross(p1 - p0, p2 - p0));
        // calculate vertex normal
        if(vert_norm)
        {
            phi0 = std::acos(glm::dot(p1 - p0, p2 - p0) / (glm::length(p1 - p0) * glm::length(p2 - p0)));
            phi1 = std::acos(glm::dot(p0 - p1, p2 - p1) / (glm::length(p0 - p1) * glm::length(p2 - p1)));
            phi2 = pi - phi0 - phi1;

            m_normals[triangle.x] += phi0 * norm;
            m_normals[triangle.y] += phi1 * norm;
            m_normals[triangle.z] += phi2 * norm;
        }
    }
    // re-normalize all normals
    for(std::vector<glm::vec3>::iterator iter = m_normals.begin(); iter != m_normals.end(); ++iter)
    {
        *iter = glm::normalize(*iter);
        m_colors.push_back(mat_color);
        //m_colors.push_back(*iter);
    }
	
	int num=m_indices.size()/3;
	vector<unsigned int> Triangles;
	for(int i=0;i<num;i++) Triangles.push_back(i);
	tree.setPosition(this->pos);
	tree.createTree(0,&m_positions,&m_indices,Triangles);
}