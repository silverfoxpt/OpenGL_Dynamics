#ifndef _TRIANGLE_DRAW_H_
#define _TRIANGLE_DRAW_H_

#include <iostream>
#include <glad/glad.h> 
#include <GLFW/glfw3.h>

class TriangleDraw {
    public:
        unsigned int VBO, VAO;

        TriangleDraw() {}
        TriangleDraw(const float* vertices, const unsigned int &size, unsigned int drawMode = GL_STATIC_DRAW);
        void Draw(unsigned int numVertices);
};

#endif