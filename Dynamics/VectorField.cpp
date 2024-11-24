#include "VectorField.h"

VectorField::VectorField(int rows, int cols, const glm::vec2& initialVec) {
    this->rows = rows; this->cols = cols;

    for (int i = 0; i < rows; i++) {
        this->vectorField.push_back(std::vector<glm::vec2>());

        for (int j = 0; j < cols; j++) {
            this->vectorField[i].push_back(initialVec);
        }
    }
}

glm::vec2 VectorField::GetVector(int rowIdx, int colIdx) {
    return this->vectorField[rowIdx][colIdx];
}

glm::vec2 VectorField::GetUnitVector(int rowIdx, int colIdx) {
    return glm::normalize(this->vectorField[rowIdx][colIdx]);
}

float VectorField::GetStrengthVector(int rowIdx, int colIdx) {
    return glm::length(this->vectorField[rowIdx][colIdx]);
}

std::vector<float> VectorField::GeneratePositionField(float startX, float startY, float squareLength) {
    std::vector<float> positionData;  // This will hold the start and end points of each vector.

    // Loop over each grid square.
    for (int i = 0; i < this->rows; i++) {
        for (int j = 0; j < this->cols; j++) {
            // Get the vector (direction and magnitude) at the current grid cell.
            glm::vec2 vector = this->vectorField[i][j];

            // Calculate the start position of the vector (bottom-left corner of the square).
            float startXPos = startX + j * squareLength;
            float startYPos = startY + i * squareLength;

            // Calculate the end position of the vector using the vector's direction and magnitude.
            glm::vec2 endVector = glm::normalize(vector) * (glm::length(vector) * squareLength / 2);

            // The end position is the start position plus the vector's direction scaled by its magnitude.
            float endXPos = startXPos + endVector.x;
            float endYPos = startYPos + endVector.y;

            // Add the start and end points of this vector to the positionData array.
            positionData.push_back(startXPos);
            positionData.push_back(startYPos);
            positionData.push_back(endXPos);
            positionData.push_back(endYPos);
        }
    }

    return positionData;
}