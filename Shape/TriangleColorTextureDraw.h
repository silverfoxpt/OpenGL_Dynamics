#ifndef _TRIANGLE_COLOR_TEXTURE_DRAW_H_
#define _TRIANGLE_COLOR_TEXTURE_DRAW_H_

#include <iostream>
#include <glad/glad.h> 
#include <GLFW/glfw3.h>

class TriangleColorTextureDraw {
    public:
        unsigned int VBO, VAO;

        TriangleColorTextureDraw() {}
        TriangleColorTextureDraw(const float* vertices, const unsigned int &size, unsigned int drawMode = GL_STATIC_DRAW);
        void Draw(unsigned int numVertices);
};

#endif