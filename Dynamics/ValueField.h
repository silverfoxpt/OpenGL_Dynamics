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
        std::vector<std::vector<float>> valueField;
        int rows, cols;
        float maxValue = 1.0f;
        glm::vec3 defaultMaxColor = glm::vec3(1.0f, 1.0f, 1.0f); //White

        std::vector<float> colorData;

        ValueField() { }
        ValueField(int rows, int cols, float initialVal);

        float GetValue(int rowIdx, int colIdx);
        void SetValue(int rowIdx, int colIdx, float newVal);
        
        std::vector<float> GenerateColorField();

        // Swap function
        void SwapWith(ValueField& other) {
            std::swap(valueField, other.valueField);
            std::swap(rows, other.rows);
            std::swap(cols, other.cols);
            std::swap(maxValue, other.maxValue);
            std::swap(defaultMaxColor, other.defaultMaxColor);
            std::swap(colorData, other.colorData);
        }

        std::vector<float> GenerateVectorPositionField(
             ValueField& vField,       // staggered v-component
            float startX, float startY,
            float squareLength,
            float maxStrength = 60.0f) ;
    
        // ---------------------------------------------------------------------
        // Per-vertex colours on the same red➞green ramp you had before.
        // ---------------------------------------------------------------------
        std::vector<float> GenerateVectorColorField(
             ValueField& vField,
            float maxStrength           = 60.0f) ;

        std::vector<float> toFlat() const
        {
            std::vector<float> flat;
            flat.reserve(rows * cols);
            for (const auto& row : valueField)
                flat.insert(flat.end(), row.begin(), row.end());
            return flat;          // NRVO / move-elided, zero extra copy outside
        }
};

#endif