#include "shader.h"


shader::shader(){}

void shader::InitShader(const char* vert_path, const char* frag_path)
{
    // create shaders and shader program
    m_vert_handle = glCreateShader(GL_VERTEX_SHADER);
    m_frag_handle = glCreateShader(GL_FRAGMENT_SHADER);
    m_shaderprog_handle = glCreateProgram();

    // load shader source from file
    const char* vert_source = textFileRead(vert_path);
    const char* frag_source = textFileRead(frag_path);
    glShaderSource(m_vert_handle, 1, &vert_source, NULL);
    glShaderSource(m_frag_handle, 1, &frag_source, NULL);
    glCompileShader(m_vert_handle);
    glCompileShader(m_frag_handle);

    // compile shader source
    GLint compiled;
    glGetShaderiv(m_vert_handle, GL_COMPILE_STATUS, &compiled);
    if(!compiled)
        printShaderInfoLog(m_vert_handle);
    glGetShaderiv(m_frag_handle, GL_COMPILE_STATUS, &compiled);
    if(!compiled)
        printShaderInfoLog(m_frag_handle);

    // bind attribute locations for the shaders
    // 0 for position, 1 for color, 2 for normal.
    glBindAttribLocation(m_shaderprog_handle, 0, "v_position");
    glBindAttribLocation(m_shaderprog_handle, 1, "v_color");
    glBindAttribLocation(m_shaderprog_handle, 2, "v_normal");
    glBindAttribLocation(m_shaderprog_handle, 3, "v_texcoord");

    // attach shader to the shader program
    glAttachShader(m_shaderprog_handle, m_vert_handle);
    glAttachShader(m_shaderprog_handle, m_frag_handle);
    glLinkProgram(m_shaderprog_handle);
    GLint linked;
    glGetProgramiv(m_shaderprog_handle, GL_LINK_STATUS, &linked);
    if(!linked)
        printLinkInfoLog(m_shaderprog_handle);

    // query uniform locations from openGL.
	
    m_vbo_handle.m_uniform_modelview = glGetUniformLocation(m_shaderprog_handle, "u_modelviewMatrix");
    m_vbo_handle.m_uniform_projection = glGetUniformLocation(m_shaderprog_handle, "u_projMatrix");
    m_vbo_handle.m_uniform_transformation = glGetUniformLocation(m_shaderprog_handle, "u_transformMatrix");
    m_vbo_handle.m_uniform_enable_texture = glGetUniformLocation(m_shaderprog_handle, "u_choose_tex");
    m_vbo_handle.m_uniform_texture_sampler = glGetUniformLocation(m_shaderprog_handle, "u_sampler1");
	m_vbo_handle.m_camera_position = glGetUniformLocation(m_shaderprog_handle, "u_camera_position");

    // activate the shader program.
    glUseProgram(m_shaderprog_handle);
}

char* shader::textFileRead(const char* fileName) 
{
    char* text;

    if (fileName != NULL) {
        FILE *file = fopen(fileName, "rt");

        if (file != NULL) {
            fseek(file, 0, SEEK_END);
            int count = ftell(file);
            rewind(file);

            if (count > 0) {
                text = (char*)malloc(sizeof(char) * (count + 1));
                count = fread(text, sizeof(char), count, file);
                text[count] = '\0';    //cap off the string with a terminal symbol, fixed by Cory
            }
            fclose(file);
        }
    }
    return text;
}

void shader::printShaderInfoLog(int shader)
{
    int infoLogLen = 0;
    int charsWritten = 0;
    GLchar *infoLog;

    glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &infoLogLen);

    // should additionally check for OpenGL errors here

    if (infoLogLen > 0)
    {
        infoLog = new GLchar[infoLogLen];
        // error check for fail to allocate memory omitted
        glGetShaderInfoLog(shader,infoLogLen, &charsWritten, infoLog);
        std::cout << "InfoLog:" << std::endl << infoLog << std::endl;
        delete [] infoLog;
    }

    // should additionally check for OpenGL errors here
}

void shader::printLinkInfoLog(int prog) 
{
    int infoLogLen = 0;
    int charsWritten = 0;
    GLchar *infoLog;

    glGetProgramiv(prog, GL_INFO_LOG_LENGTH, &infoLogLen);

    // should additionally check for OpenGL errors here

    if (infoLogLen > 0)
    {
        infoLog = new GLchar[infoLogLen];
        // error check for fail to allocate memory omitted
        glGetProgramInfoLog(prog,infoLogLen, &charsWritten, infoLog);
        std::cout << "InfoLog:" << std::endl << infoLog << std::endl;
        delete [] infoLog;
    }
}

void shader::ActivateShaderprog()
{
    GLint current_prog;
    glGetIntegerv(GL_CURRENT_PROGRAM, &current_prog);
    if(current_prog != (GLint)m_shaderprog_handle)
        glUseProgram(m_shaderprog_handle);
}

void shader::DeactivateShaderprog()
{
    GLint current_prog;
    glGetIntegerv(GL_CURRENT_PROGRAM, &current_prog);
    if(current_prog == (GLint)m_shaderprog_handle)
        glUseProgram(0);
}

void shader::CleanupShader()
{
    glDetachShader(m_shaderprog_handle, m_vert_handle);
    glDetachShader(m_shaderprog_handle, m_frag_handle);
    glDeleteShader(m_vert_handle);
    glDeleteShader(m_frag_handle);
    glDeleteProgram(m_shaderprog_handle);
}

void shader::SetCameraProjection(glm::mat4 projection)
{
    glMatrixMode(GL_PROJECTION);
    glLoadMatrixf(&projection[0][0]);
    
    ActivateShaderprog();
    glUniformMatrix4fv(m_vbo_handle.m_uniform_projection, 1, false, &projection[0][0]);
}

void shader::SetCameraModelview(glm::mat4 modelview)
{
    glMatrixMode(GL_MODELVIEW);
    glLoadMatrixf(&modelview[0][0]);
    
    ActivateShaderprog();
    glUniformMatrix4fv(m_vbo_handle.m_uniform_modelview, 1, false, &modelview[0][0]);
}