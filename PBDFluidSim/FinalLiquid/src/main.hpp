#ifndef _MAIN_
#define _MAIN_

#include <cstdlib>
#include "globalHeader.h"
#include "constriant.h"

#include "scene.h"
#include "simulation.h"
#include "container.h"
#include "shader.h"
#include "io.h"
#include "kernel.h"

glm::mat4x4 getCameraProjectionView();
glm::mat4x4 getCameraModelView();

#endif