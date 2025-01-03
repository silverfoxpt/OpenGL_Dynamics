#ifndef _LINES_COLOR_DRAW_H_
#define _LINES_COLOR_DRAW_H_

#include <iostream>
#include <glad/glad.h> 
#include <GLFW/glfw3.h>
#include <vector>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

class LinesColorDraw {
public:
    unsigned int VBO_position, VBO_color, VAO;
    std::vector<float> colors, positions;

    LinesColorDraw();
    LinesColorDraw(const std::vector<float>& position, const std::vector<float>& color);

    void UpdatePositions(const std::vector<float>& newPositions);
    void UpdateColors(const std::vector<float>& newColors);
    void Draw(float lineWidth);
};

#endif