#include "Particle.h"

Particle::Particle(){
	this->pos=glm::vec3(0,0,0);
}

Particle::Particle(float radius){
	this->radius=radius;
	this->pos=glm::vec3(0,0,0);
	this->vel=glm::vec3(0,0,0);
	this->acc=glm::vec3(0,0,0);
	this->dist=glm::vec3(0,0,0);
}

Particle::Particle(float radius,glm::vec3 pos){
	this->radius=radius;
	this->pos=pos;
	this->vel=glm::vec3(0,0,0);
	this->acc=glm::vec3(0,0,0);
	this->dist=glm::vec3(0,0,0);
	this->init();
}

void Particle::draw(const VBO& vbos){
	
	glEnable(GL_LIGHTING);
	glPushMatrix();
	glTranslated(pos[0],pos[1],pos[2]);
	glutSolidSphere(radius, 10, 10);
	glPopMatrix();
	glDisable(GL_LIGHTING);
	
	/*
	glBindBuffer(GL_ARRAY_BUFFER, vbos.m_vbo);
    glBufferData(GL_ARRAY_BUFFER, 3 * m_positions.size() * sizeof(float), &m_positions[0], GL_DYNAMIC_DRAW);
    // color
    glBindBuffer(GL_ARRAY_BUFFER, vbos.m_cbo);
    glBufferData(GL_ARRAY_BUFFER, 3 * m_colors.size() * sizeof(float), &m_colors[0], GL_DYNAMIC_DRAW);

    // normal
    glBindBuffer(GL_ARRAY_BUFFER, vbos.m_nbo);
    glBufferData(GL_ARRAY_BUFFER, 3 * m_normals.size() * sizeof(float), &m_normals[0], GL_DYNAMIC_DRAW);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbos.m_ibo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, m_indices.size() * sizeof(unsigned short), &m_indices[0], GL_DYNAMIC_DRAW);
	
    glEnableVertexAttribArray(0);
    glEnableVertexAttribArray(1);
    glEnableVertexAttribArray(2);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbos.m_ibo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, m_indices.size() * sizeof(unsigned short), &m_indices[0], GL_DYNAMIC_DRAW);
	
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
	*/
}

/*
brought from tiantianLiu's sphere visualization
*/
void Particle::init(){
    m_positions.clear();
    m_colors.clear();
    m_normals.clear();
    m_indices.clear();

    glm::vec3 mat_color(0.3,0.5,0.8);
    unsigned int slice = 12, stack = 5;

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
/*
brought from tiantianLiu's sphere visualization
*/

float Particle::getRadius(){
	return radius;
}

float Particle::getRest(){
	return rest;
}