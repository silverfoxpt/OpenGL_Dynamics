#include "TriangleDraw.h"

// Original, as a reminder if needed
// // Gen VAO - All subsequent calls on setup to VBO will be saved to VAO - reusing available
// glGenVertexArrays(1, &triangleVAO);
// glBindVertexArray(triangleVAO);

// // Gen VBO
// glGenBuffers(1, &triangleVBO);
// glBindBuffer(GL_ARRAY_BUFFER, triangleVBO);

// // Setup VBO
// glBufferData(GL_ARRAY_BUFFER, sizeof(triangleVertex), triangleVertex, GL_STATIC_DRAW);
// glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*) 0);
// glEnableVertexAttribArray(0);

// // Reset bind
// glBindBuffer(GL_ARRAY_BUFFER, 0);
// glBindVertexArray(0);

TriangleDraw::TriangleDraw(const float* vertices, const unsigned int &size, unsigned int drawMode) {
    glGenVertexArrays(1, &this->VAO);
    glBindVertexArray(this->VAO);

    glGenBuffers(1, &this->VBO);
    glBindBuffer(GL_ARRAY_BUFFER, this->VBO);

    // Load vertex data into VBO
    glBufferData(GL_ARRAY_BUFFER, size, vertices, drawMode);

    // Set vertex attributes
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*) 0);
    glEnableVertexAttribArray(0);

    // Unbind VBO and VAO
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
}

void TriangleDraw::Draw(unsigned int numVertices) {
    glBindVertexArray(this->VAO);
    glDrawArrays(GL_TRIANGLES, 0, numVertices);
    glBindVertexArray(0);
}