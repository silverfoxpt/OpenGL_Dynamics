#include "ValueField.h"

ValueField::ValueField(int rows, int cols, float initialVal) {
    this->rows = rows; this->cols = cols;
    if (glm::length(initialVal) > this->maxValue) {
        initialVal = this->maxValue;
    }

    for (int i = 0; i < rows; i++) {
        this->valueField.push_back(std::vector<float>());

        for (int j = 0; j < cols; j++) {
            this->valueField[i].push_back(initialVal);
        }
    }
}

float ValueField::GetValue(int rowIdx, int colIdx) {
    return this->valueField[rowIdx][colIdx];
}

void ValueField::SetValue(int rowIdx, int colIdx, float newVal) {
    this->valueField[rowIdx][colIdx] = newVal;
}

std::vector<float> ValueField::GenerateColorField() {
    std::vector<float> colorData;  // This will hold the color of each arrow.
    colorData.reserve(this->rows * this->cols * 6);

    // Loop over each grid square.
    for (int i = 0; i < this->rows; i++) {
        for (int j = 0; j < this->cols; j++) {
            // Get the value at the current grid cell.
            float value = GetValue(i, j);
            float strength = value / this->maxValue;

            // Add color to color array
            for (int k = 0; k <  6; k++) {
                colorData.push_back(defaultMaxColor.r * strength);
                colorData.push_back(defaultMaxColor.g * strength);
                colorData.push_back(defaultMaxColor.b * strength);
            }
            
        }
    }

    return colorData;
}