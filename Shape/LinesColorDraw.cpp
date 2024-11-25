#include "LinesColorDraw.h"

LinesColorDraw::LinesColorDraw() {}

LinesColorDraw::LinesColorDraw(const std::vector<float>& verticesPos, const std::vector<float>& verticesColor) {
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

void LinesColorDraw::Draw(float lineWidth) {
    glLineWidth(lineWidth);
    glBindVertexArray(VAO);
    glDrawArrays(GL_LINES, 0, 2); // Draw the arrow as a line
    glBindVertexArray(0);
}
