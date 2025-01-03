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
    colorData.resize(this->rows * this->cols * 6 * 3);
}

float ValueField::GetValue(int rowIdx, int colIdx) {
    return this->valueField[rowIdx][colIdx];
}

void ValueField::SetValue(int rowIdx, int colIdx, float newVal) {
    this->valueField[rowIdx][colIdx] = newVal;
}

std::vector<float> ValueField::GenerateColorField() {
    int cnt = 0;

    // Loop over each grid square
    for (int i = 0; i < this->rows; i++) {
        for (int j = 0; j < this->cols; j++) {
            // Get the value at the current grid cell
            float value = GetValue(i, j);
            float strength = value / this->maxValue;

            // Precompute scaled color values
            float r = defaultMaxColor.r * strength;
            float g = defaultMaxColor.g * strength;
            float b = defaultMaxColor.b * strength;

            // Add color to color array
            for (int k = 0; k < 6; k++) {
                colorData[cnt++] = r;
                colorData[cnt++] = g; 
                colorData[cnt++] = b; 
            }
            
        }
    }
    return colorData;
}