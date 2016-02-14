# PBDFluid
## This is the project for Penn CIS563(Physical-based Animation)
## Team Member: Siqi Huang, Hansi Liu, and Weijian Zhou

# Real-time Demo(20*25*20 Particles Collision with objects)
## CPU(Inter(R) Core(TM) i7-4870HQ 2.5GHZ)
[![ScreenShot](pic/screenshot.png)](https://youtu.be/hTkOVZNK4sQ)

## GPU(NVIDIA GeForce GTX 970M)
[![ScreenShot](pic/screenshot.png)](https://youtu.be/-JkCEV8yWuM)

Those two videos basically rendered the same scene, but obviously after using GPU to accelerate, the render is more than 10 times faster(of course the actual speed depend on your own hardware).

# Part 0 Architecture and Control(Siqi Huang)
## Architecture
In this project, we have used some other's code as referrence, which it is necessary to acknowledge here:
1 The basic openGl view is gluv and glew. For gluv part, we have using some of Zonhan Xu's code as referrence.
2 The basic structure of the code is referred to Tiantian Liu's Deformable_body_Sim
3 The gl buffer part, vertex and fragment shader is mostly from Tiantian Liu's project

## Control
### Real-time input
use mouse to rotate the scene, the gravity change's with the view.
use blankspace to move to the next frame

### Input
use CPU_RENDER or GPU_RENDER to specify render mode<br>
use CONTAINER_CUBE or CONTAINER_SPHERE to specify container shape(container always at origin) and x y z(or r) to specify dimension(or radius)<br>
use EYE,LOOKAT,UP to setup camera<br>
use PARTICLE x y z to specify particle number<br>
use MINPOS to specify start position of particles<br>
use DIS to specify particle distance<br>
