#include "VectorField.h"

VectorField::VectorField(int rows, int cols, glm::vec2 initialVec) {
    this->rows = rows; this->cols = cols;
    if (glm::length(initialVec) > this->maxStrength) {
        initialVec = glm::normalize(initialVec) * maxStrength;
    }

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

void VectorField::SetVector(int rowIdx, int colIdx, glm::vec2 newValue) {
    this->vectorField[rowIdx][colIdx] = newValue;
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
            float startXPos = startX + j * squareLength + squareLength / 2;
            float startYPos = startY - i * squareLength - squareLength / 2;

            // Calculate the end position of the vector
            glm::vec2 endVector = glm::normalize(vector) * (squareLength / 2);

            // The end position is the start position plus the vector's direction scaled by its magnitude.
            float endXPos = startXPos + endVector.x;
            float endYPos = startYPos + endVector.y;

            // Add the start and end points of this vector to the positionData array.
            positionData.push_back(startXPos);
            positionData.push_back(startYPos);
            positionData.push_back(0);

            positionData.push_back(endXPos);
            positionData.push_back(endYPos);
            positionData.push_back(0);
        }
    }

    return positionData;
}

std::vector<float> VectorField::GenerateColorField() {
    std::vector<float> colorData;  // This will hold the start and end points of each vector.

    // Loop over each grid square.
    for (int i = 0; i < this->rows; i++) {
        for (int j = 0; j < this->cols; j++) {
            // Get the vector (direction and magnitude) at the current grid cell.
            glm::vec2 vector = this->vectorField[i][j];
            float strength = glm::length(vector) / this->maxStrength;

            // Add the start and end points of this vector to the positionData array.
            colorData.push_back(defaultStartColor.r * strength);
            colorData.push_back(defaultStartColor.g * strength);
            colorData.push_back(defaultStartColor.b * strength);

            colorData.push_back(defaultEndColor.r * strength);
            colorData.push_back(defaultEndColor.g * strength);
            colorData.push_back(defaultEndColor.b * strength);
        }
    }

    return colorData;
}

void VectorField::TestUpdate() {
    float angleDegrees = -0.1f;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            glm::vec2 curVec = GetVector(i, j);

            // Convert degrees to radians
            float angleRadians = glm::radians(angleDegrees);

            // Calculate the rotated coordinates
            float x = curVec.x * glm::cos(angleRadians) - curVec.y * glm::sin(angleRadians);
            float y = curVec.x * glm::sin(angleRadians) + curVec.y * glm::cos(angleRadians);

            SetVector(i, j, glm::vec2(x, y));
        }
    }
}
