#ifndef _SHADER_
#define _SHADER_

#include "globalHeader.h"
#include "openglHeader.h"

/*
This class referred to Tiantian Liu's deformable body sim code
*/

class shader{
public:
	shader();
	void InitShader(const char* vert_path, const char* frag_path);
	char* shader::textFileRead(const char* fileName);
	void shader::printShaderInfoLog(int shader);
	void shader::printLinkInfoLog(int prog);
	void shader::ActivateShaderprog();
	void shader::DeactivateShaderprog();
	void shader::SetCameraProjection(glm::mat4 projection);
	void shader::SetCameraModelview(glm::mat4 modelview);
	void shader::CleanupShader();

	VBO m_vbo_handle;
    GLuint m_vert_handle, m_frag_handle, m_shaderprog_handle;
    GLuint m_texture;

public: // inlines
    inline VBO& getVBO() {return m_vbo_handle;}
};

#endif