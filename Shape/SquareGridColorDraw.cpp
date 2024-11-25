#include "SquareGridColorDraw.h"

std::vector<float> SquareGridColorDraw::GenerateSampleGrid(int rows, int cols, float squareSize) {
    std::vector<float> vertices;
    for (int row = 0; row < rows; ++row) {
        for (int col = 0; col < cols; ++col) {
            float x = col * squareSize;
            float y = -row * squareSize;

            // Define the two triangles for this square
            float squareVertices[] = {
                // Triangle 1
                x, y, 0.0f,
                x + squareSize, y, 0.0f,
                x, y - squareSize, 0.0f,

                // Triangle 2
                x + squareSize, y, 0.0f,
                x + squareSize, y - squareSize, 0.0f,
                x, y - squareSize, 0.0f
            };
            vertices.insert(vertices.end(), std::begin(squareVertices), std::end(squareVertices));
        }
    }

    return vertices;
}

std::vector<float> SquareGridColorDraw::GenerateSampleGridColors(int rows, int cols) {
    std::vector<float> vertices;
    for (int row = 0; row < rows; ++row) {
        for (int col = 0; col < cols; ++col) {
            // Define the two triangles' color for this square
            float squareVertices[] = {
                // Triangle 1
                0.617f, 0.617f, 0.617f,
                0.617f, 0.617f, 0.617f,
                0.617f, 0.617f, 0.617f,

                // Triangle 2
                0.617f, 0.617f, 0.617f,
                0.617f, 0.617f, 0.617f,
                0.617f, 0.617f, 0.617f
            };
            vertices.insert(vertices.end(), std::begin(squareVertices), std::end(squareVertices));
        }
    }

    return vertices;
}

SquareGridColorDraw::SquareGridColorDraw(
    int rows, int cols, float size,
    const std::vector<float>& verticesPos, const std::vector<float>& verticesColor, unsigned int drawMode) 
{
    this->rows = rows; this->cols = cols; this->size = size;

    glGenVertexArrays(1, &this->VAO);
    glBindVertexArray(this->VAO);

    // Gen position VBO
    glGenBuffers(1, &this->VBO_position);
    glBindBuffer(GL_ARRAY_BUFFER, this->VBO_position);

    // Load vertex data into VBO
    glBufferData(GL_ARRAY_BUFFER, verticesPos.size() * sizeof(float), verticesPos.data(), drawMode);

    // Set vertex attributes
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*) 0);
    glEnableVertexAttribArray(0);

    // Gen color VBO
    glGenBuffers(1, &this->VBO_color);
    glBindBuffer(GL_ARRAY_BUFFER, this->VBO_color);

    // Load vertex data into VBO
    glBufferData(GL_ARRAY_BUFFER, verticesColor.size() * sizeof(float), verticesColor.data(), drawMode);

    // Set vertex attributes
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, (void*) 0);
    glEnableVertexAttribArray(1);    

    // Unbind VBO and VAO
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
}

void SquareGridColorDraw::UpdateColors(const std::vector<float>& colors) {
    glBindBuffer(GL_ARRAY_BUFFER, VBO_color);
    glBufferSubData(GL_ARRAY_BUFFER, 0, colors.size() * sizeof(float), colors.data());
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void SquareGridColorDraw::Draw(unsigned int numVertices) {
    glBindVertexArray(this->VAO);
    glDrawArrays(GL_TRIANGLES, 0, numVertices);
    glBindVertexArray(0);
}