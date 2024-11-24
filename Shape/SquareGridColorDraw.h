#ifndef _SQUARE_GRID_COLOR_DRAW_H_
#define _SQUARE_GRID_COLOR_DRAW_H_

#include <iostream>
#include <glad/glad.h> 
#include <GLFW/glfw3.h>
#include <vector>

class SquareGridColorDraw {
    public:
        unsigned int VBO_position, VBO_color, VAO;
        int rows, cols;
        float size;

        SquareGridColorDraw() {}
        SquareGridColorDraw(
            int rows, int cols, float size,
            const std::vector<float>& verticesPos, const std::vector<float>& verticesColor, unsigned int drawMode = GL_STATIC_DRAW
        );
        void Draw(unsigned int numVertices);
        void UpdateColors(const std::vector<float>& colors);
        
        static std::vector<float> GenerateSampleGrid(int rows, int cols, float squareSize);
        static std::vector<float> GenerateSampleGridColors(int rows, int cols);
};

#endif