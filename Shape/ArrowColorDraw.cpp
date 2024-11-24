#include "ArrowColorDraw.h"

ArrowColorDraw::ArrowColorDraw() {}

ArrowColorDraw::ArrowColorDraw(
    const glm::vec3& start, const glm::vec3& end, const glm::vec3& initialColor
) {
    this->startPosition = start; this->endPosition = end; this->color = initialColor;

    // Initialize vertex data
    std::vector<float> verticesPos = {
        start.x, start.y, start.z,
        end.x, end.y, end.z
    };
    std::vector<float> verticesColor = {
        color.r, color.g, color.b,
        color.r, color.g, color.b
    };

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

void ArrowColorDraw::UpdatePositions(const glm::vec3& newStart, const glm::vec3& newEnd) {
    startPosition = newStart;
    endPosition = newEnd;

    std::vector<float> verticesPos = {
        newStart.x, newStart.y, newStart.z,
        newEnd.x, newEnd.y, newEnd.z
    };

    glBindBuffer(GL_ARRAY_BUFFER, VBO_position);
    glBufferSubData(GL_ARRAY_BUFFER, 0, verticesPos.size() * sizeof(float), verticesPos.data());
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void ArrowColorDraw::UpdateColors(const glm::vec3& newColor) {
    color = newColor;

    std::vector<float> verticesColor = {
        newColor.r, newColor.g, newColor.b,
        newColor.r, newColor.g, newColor.b
    };

    glBindBuffer(GL_ARRAY_BUFFER, VBO_color);
    glBufferSubData(GL_ARRAY_BUFFER, 0, verticesColor.size() * sizeof(float), verticesColor.data());
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void ArrowColorDraw::Draw(float lineWidth) {
    glLineWidth(lineWidth);
    glBindVertexArray(VAO);
    glDrawArrays(GL_LINES, 0, 2); // Draw the arrow as a line
    glBindVertexArray(0);
}
