#ifndef _TRIANGLE_COLOR_DRAW_H_
#define _TRIANGLE_COLOR_DRAW_H_

#include <iostream>
#include <glad/glad.h> 
#include <GLFW/glfw3.h>

class TriangleColorDraw {
    public:
        unsigned int VBO, VAO;

        TriangleColorDraw() {}
        TriangleColorDraw(const float* vertices, const unsigned int &size, unsigned int drawMode = GL_STATIC_DRAW);
        void Draw(unsigned int numVertices);
};

#endif