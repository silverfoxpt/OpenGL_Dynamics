#include "LinesColorDraw.h"

LinesColorDraw::LinesColorDraw() {}

LinesColorDraw::LinesColorDraw(const std::vector<float>& verticesPos, const std::vector<float>& verticesColor) {
    this->colors = verticesColor; this->positions = verticesPos;
    
    // Generate buffers and array object
    glGenVertexArrays(1, &VAO);
    glBindVertexArray(VAO);

    glGenBuffers(1, &VBO_position);

    // Position buffer
    glBindBuffer(GL_ARRAY_BUFFER, VBO_position);
    glBufferData(GL_ARRAY_BUFFER, verticesPos.size() * sizeof(float), verticesPos.data(), GL_DYNAMIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
    glEnableVertexAttribArray(0);

    glGenBuffers(1, &VBO_color);

    // Color buffer
    glBindBuffer(GL_ARRAY_BUFFER, VBO_color);
    glBufferData(GL_ARRAY_BUFFER, verticesColor.size() * sizeof(float), verticesColor.data(), GL_DYNAMIC_DRAW);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
    glEnableVertexAttribArray(1);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
}

void LinesColorDraw::UpdatePositions(const std::vector<float>& newPositions) {
    this->positions = newPositions;
    glBindBuffer(GL_ARRAY_BUFFER, VBO_position);
    glBufferSubData(GL_ARRAY_BUFFER, 0, newPositions.size() * sizeof(float), newPositions.data());
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void LinesColorDraw::UpdateColors(const std::vector<float>& newColors) {
    this->colors = newColors;
    glBindBuffer(GL_ARRAY_BUFFER, VBO_color);
    glBufferSubData(GL_ARRAY_BUFFER, 0, newColors.size() * sizeof(float), newColors.data());
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void LinesColorDraw::Draw(float lineWidth) {
    glLineWidth(lineWidth);

    glBindVertexArray(VAO);
    glDrawArrays(GL_LINES, 0, this->positions.size() / 3); // Draw the arrow as a line
    glBindVertexArray(0);
}
