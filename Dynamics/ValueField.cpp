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

#include "ValueField.h"
#include <algorithm>   // std::clamp
#include <cmath>       // std::sqrt

// ------------------------------------------------------------
// Arrow vertices (cell centre → tip)
// ------------------------------------------------------------
std::vector<float> ValueField::GenerateVectorPositionField(
         ValueField& vField,
        float startX, float startY,
        float squareLength,
        float maxStrength) 
{
    // sanity check – both fields must be same size
    if (rows != vField.rows || cols != vField.cols)
        throw std::runtime_error("u and v fields dimension mismatch");

    std::vector<float> data;
    data.reserve(rows * cols * 6);

    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            glm::vec2 vel{ GetValue(i, j), vField.GetValue(i, j) };
            float strength = glm::length(vel) / maxStrength;

            float cx = startX + j * squareLength + squareLength * 0.5f;
            float cy = startY - i * squareLength - squareLength * 0.5f;

            glm::vec2 tip = glm::normalize(vel) *
                            (squareLength * 0.5f) * strength;

            data.insert(data.end(),
            {   cx,           cy,           0.0f,      // tail
                cx + tip.x,   cy + tip.y,   0.0f });   // head
        }
    }
    return data;
}

// ------------------------------------------------------------
// Arrow colours (duplicated for head & tail)
// ------------------------------------------------------------
std::vector<float> ValueField::GenerateVectorColorField(
         ValueField& vField,
        float maxStrength) 
{
    glm::vec3 startColor = {1.0f, 0.0f, 0.0f};
    glm::vec3 endColor = {0.0f, 1.0f, 0.0f};

    if (rows != vField.rows || cols != vField.cols)
        throw std::runtime_error("u and v fields dimension mismatch");

    std::vector<float> colors;
    colors.reserve(rows * cols * 6);

    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            glm::vec2 vel{ GetValue(i, j), vField.GetValue(i, j) };
            float t = std::clamp(glm::length(vel) / maxStrength, 0.0f, 1.0f);

            glm::vec3 c = (1.0f - t) * startColor + t * endColor;

            // duplicate for both vertices of the arrow
            for (int vtx = 0; vtx < 2; ++vtx)
            {
                colors.push_back(c.r);
                colors.push_back(c.g);
                colors.push_back(c.b);
            }
        }
    }
    return colors;
}
