
# 3D Computer Graphics assignment ([Lab2](https://hackmd.io/5O2FIpo7RuCUEnjf0qvQiA?view))
### Work requirements:
1. There are many test data in the forder and many input commands in each test data. Please read the files and complete the following instructions.
1. You can only use the function of drawing **points** in OpenGL (other functions such as line drawing cannot be used.
1. Eigen or other outer libaries cannot be used.

## Environment
`Linux distribution`: `Ubuntu 22.04`<br>
About how to build the enviroment on linux: [here](https://hackmd.io/3xPNjv6kRh2Ml6Ll7Nlw4A).

## Operating instructions
* ### Comment:
    The instruction is `# ......`<br>
    If a `#` appears, skip the line. (i.e. You don't need to execute the line which # appears.)  
    
* ### Reset:
    The instruction is `reset`<br>
    Reset the transformation matrix.
    
* ### Translate:
    The instruction is `translate x y`<br>
    `x` means moving along x-axis by `x` unit.<br>
    `y` means moving along y-axis by `y` unit.

* ### Scale:
    The instruction is `scale x y`<br>
    `x` means scaling along x-axis by `x` unit.<br>
    `y` means scaling along y-axis by `y` unit.
    
* ### Rotate:
    The instruction is `rotate θ`<br>
    `θ` means rotating along z-axis by θ degree.

* ### Clear data:
    The instruction is `clearData`<br>
    Clear the objects that you create (square and triangle).

* ### Clear screen:
    The instruction is `clearScreen`<br>
    Clear the screen.

* ### View:
    The instruction is `view wxl wxr wyb wyt vxl vxr vyb vyt`<br>
    `wxl wxr wyb wyt` means the position of object before mapping.<br>
    `vxl vxr vyb vyt` means the position of object after mapping.<br>
    
    In math words <br>

        f(wxl,wyb) = (vxl,vyb)
        f(wxr,wyt) = (vxr,vyt)
where it is a linear transformation.<br>
The objects outside the boundary need to **clipped**.<br>
Every time you execute `view`, the created objects must be displayed. Please draw the boundary of objects with points and lines.(Use the code in lab1)
Please add `system("pause");` when you display anything.

* ### Square:
    The instruction is `square`<br>
    Create a square with its vertex `(-1,-1)` `(1,-1)` `(1,1)` `(-1,1)`.

* ### Triangle:
    The instruction is `triangle`<br>
    Create a triangle with its vertex `(0,1)` `(-1,-1)` `(1,-1)`.

* ### End:
    The instruction is `end`<br>
    Destory the window.

* ### Bonus function:
    Each object automatically has a different color and be filled.

## Libraries
* ### C++
```cpp
#include <GL/glut.h>
#include <GL/gl.h>
#include <cmath> 
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <fstream>
#include <string>
```

* ### freeGLUT
```cpp
glutInit(int *&argcp, char **argv);
glutInitWindowSize(int width, int height);
glutInitWindowPosition(int x, int y);
glutInitDisplayMode(unsigned int mode);
glutSwapBuffers(void);
glutMainLoop();
glutCreateWindow(char *name);
glutDisplayFunc(void (*func)(void));
```
```
gluOrtho2D(GLdouble left, GLdouble right, GLdouble top, GLdouble bottom); 
glClear(GLbitfield mask);
glPointSize(GLfloat size);
glColor3f(GLfloat red, GLfloat green, GLfloat blue);
glBegin(GLenum mode);
glVertex2i(GLint x, GLint y);
glEnd(void);
```
