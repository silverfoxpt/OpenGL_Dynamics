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
    glm::vec3 startPosition, endPosition, color;

    LinesColorDraw();
    LinesColorDraw(const std::vector<float>& position, const std::vector<float>& color);

    void UpdatePositions(const glm::vec3& newStart, const glm::vec3& newEnd);
    void UpdateColors(const glm::vec3& newColor);
    void Draw(float lineWidth);
};

#endif