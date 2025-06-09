#ifndef _VALUE_FIELD_H_
#define _VALUE_FIELD_H_

#include <iostream>
#include <glad/glad.h> 
#include <GLFW/glfw3.h>
#include <vector>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

class ValueField {
    public:
        std::vector<float> valueField; // Now 1D
        int rows, cols;
        float maxValue = 1.0f;
        glm::vec3 defaultMaxColor = glm::vec3(1.0f, 1.0f, 1.0f); //White

        std::vector<float> colorData;

        ValueField() : rows(0), cols(0) { }
        ValueField(int rows, int cols, float initialVal);

        inline float GetValue(int rowIdx, int colIdx) const {
            return valueField[rowIdx * cols + colIdx];
        }
        inline void SetValue(int rowIdx, int colIdx, float newVal) {
            valueField[rowIdx * cols + colIdx] = newVal;
        }

        std::vector<float> GenerateColorField();

        void SwapWith(ValueField& other) {
            std::swap(valueField, other.valueField);
            std::swap(rows, other.rows);
            std::swap(cols, other.cols);
            std::swap(maxValue, other.maxValue);
            std::swap(defaultMaxColor, other.defaultMaxColor);
            std::swap(colorData, other.colorData);
        }

        std::vector<float> GenerateVectorPositionField(
             ValueField& vField,
            float startX, float startY,
            float squareLength,
            float maxStrength = 60.0f);

        std::vector<float> GenerateVectorColorField(
             ValueField& vField,
            float maxStrength = 60.0f);

        std::vector<float> toFlat() const
        {
            return valueField;
        }
};

#endif