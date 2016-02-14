#include "main.hpp"
#include <ctime>

// GUI interaction stuff
GLVU glvu;
scene *g_scene;
simulation *g_simulation;
container *g_container;
shader *g_render;
fileInput *fileinput;
glm::vec3 gravityVector=glm::vec3(0,gravity,0);
float *cameraModel,*cameraProjection;
bool animate = false;

void init(){
	glewInit();
	g_render=new shader();
	g_render->InitShader("./shaders/vert.glsl","./shaders/frag.glsl");
	cameraModel=new float[16];
	cameraProjection=new float[16];

	switch(fileinput->getContainerType()){
	case CONTAINER::CONTAINER_CUBE:
		g_container=new container_cube(glm::vec3(0,0,0),fileinput->getCubeDim(),smooth_radius);
		break;
	case CONTAINER::CONTAINER_SPHERE:
		g_container=new container_sphere(glm::vec3(0,0,0),fileinput->getSphereRadius(),smooth_radius);
	}
	
	g_scene=new scene();
	if(!g_scene->init(fileinput,g_container))
		exit(0);
	g_simulation=new simulation(g_scene);


	/*
	GPU init
	*/
	if(fileinput->isGPURender()){
		initParticleOnGPU(fileinput->getParticleDim(),fileinput->getParticleMinPos(),fileinput->getParticleDis());
		switch(fileinput->getContainerType()){
		case CONTAINER::CONTAINER_CUBE:
			initCubeContainerOnGPU(fileinput->getCubeDim(),smooth_radius);
			break;
		case CONTAINER::CONTAINER_SPHERE:
			initSphereContainerOnGPU(fileinput->getSphereRadius(),smooth_radius);
			break;
		}
		GPUPrimitive *Gprimitive=new GPUPrimitive[fileinput->pri.size()];
		for(int i=0;i<fileinput->pri.size();++i){
			Gprimitive[i]=fileinput->pri[i]->getGPUPrimitive();
		}
		initPrimitiveOnGPU(Gprimitive,fileinput->pri.size());
	}
	/*
	GPU init
	*/
}

void displayCallback()
{
	glvu.BeginFrame();

	// clear away the previous frame
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_DEPTH_TEST);

	g_simulation->update(fileinput->isGPURender());
	
	g_scene->drawParticle(g_render->getVBO(),fileinput->isGPURender());
	g_container->draw();

	g_render->SetCameraModelview(getCameraModelView());
	g_render->SetCameraProjection(getCameraProjectionView());
	g_scene->drawPrimitive(g_render->getVBO());
	g_render->ActivateShaderprog();
	g_render->DeactivateShaderprog();
	
	glvu.EndFrame();
}


void keyboardCallback(unsigned char key, int x, int y)
{
	glutPostRedisplay();
}

void glutMouseClick(int button, int state, int x, int y)
{
	glvu.Mouse(button,state,x,y);
}

void glutMouseMotion(int x, int y)
{
	glvu.Motion(x,y);
	glvuVec3f viewVector = glvu.GetCurrentCam()->Y;
	gravityVector = gravity * glm::vec3(viewVector[0],viewVector[1],viewVector[2]);
}


void idleCallback()
{
	if (!animate) return;
	glutPostRedisplay();
}

glm::mat4x4 getCameraModelView(){
	cameraModel=glvu.GetModelviewMatrix(cameraModel);
	glm::mat4x4 m;
	for(int i=0;i<4;++i){
		for(int j=0;j<4;++j){
			m[i][j]=cameraModel[i*4+j];
		}
	}
	return m;
}

glm::mat4x4 getCameraProjectionView(){
	cameraProjection=glvu.GetProjectionMatrix(cameraProjection);
	glm::mat4x4 m;
	for(int i=0;i<4;++i){
		for(int j=0;j<4;++j){
			m[i][j]=cameraProjection[i*4+j];
		}
	}
	return m;
}

int main(int argc, char** argv)
{
	srand((unsigned)time(0));  
	fileinput=new fileInput(DEFAULT_SCENE_FILE);
	if(!fileinput->readFile()) exit(0);

	char title[40]="PBD Fluid Simulation(Rendered On CPU)";
	if(fileinput->isGPURender()){ 
		title[33]='G';
	}

	glutInit(&argc, argv);
	glvu.Init(title, GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH, 0, 0, 800, 800);
	glShadeModel(GL_SMOOTH);

	glvu.SetInertiaEnabled(0);

	// point GLUT to our callback functions
	glutDisplayFunc(displayCallback); 
	glutIdleFunc(idleCallback);
	glutKeyboardFunc(keyboardCallback);
	glutMouseFunc(glutMouseClick);
	glutMotionFunc(glutMouseMotion);

	init();

	// set background to white
	glClearColor(1.0, 1.0, 1.0, 1.0);

	// enable lights
	GLfloat ambient[] = {0.7,0.7,0.7};
	GLfloat diffuse[] = {1.0,1.0,1.0};
	GLfloat specular[] = {0.0, 0.0, 0.0};

	glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, specular);

	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	

	glvuVec3f ModelMin(-10,-10,-10), ModelMax(10,10,10);
	
	float Yfov = 45;
	float Aspect = 1;
	float Near = 0.001f;
	float Far = 10.0f;
	glvu.SetAllCams(ModelMin,ModelMax,fileinput->getCameraPos(),fileinput->getCameraLookAt(),fileinput->getCameraUp(),Yfov,Aspect,Near,Far);

	glvuVec3f center(0.0, 0.0, 0.0);
	glvu.SetWorldCenter(center);

	// Let GLUT take over
	glutMainLoop();
	getchar();
	return 0;
}
