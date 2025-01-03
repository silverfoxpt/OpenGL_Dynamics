#ifndef _VECTOR_FIELD_H_
#define _VECTOR_FIELD_H_

#include <iostream>
#include <glad/glad.h> 
#include <GLFW/glfw3.h>
#include <vector>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

class VectorField {
    public:
        std::vector<std::vector<glm::vec2>> vectorField;
        int rows, cols;
        float maxStrength = 60.0f;
        glm::vec3 defaultStartColor = glm::vec3(1.0f, 0.0f, 0.0f); //Red
        glm::vec3 defaultEndColor = glm::vec3(0.0f, 1.0f, 0.0f); //Red

        VectorField() { }
        VectorField(int rows, int cols, glm::vec2 initialVec);

        glm::vec2 GetVector(int rowIdx, int colIdx);
        void SetVector(int rowIdx, int colIdx, glm::vec2 newValue);

        glm::vec2 GetUnitVector(int rowIdx, int colIdx);
        float GetStrengthVector(int rowIdx, int colIdx);

        std::vector<float> GeneratePositionField(float startX, float startY, float length);
        std::vector<float> GenerateColorField();
        void TestUpdate();

        // Swap function
        void SwapWith(VectorField& other) {
            std::swap(vectorField, other.vectorField);
            std::swap(rows, other.rows);
            std::swap(cols, other.cols);
            std::swap(maxStrength, other.maxStrength);
            std::swap(defaultStartColor, other.defaultStartColor);
            std::swap(defaultEndColor, other.defaultEndColor);
        }

        // CopyFrom function
        void CopyFrom(const VectorField& other) {
            vectorField = other.vectorField;
            rows = other.rows;
            cols = other.cols;
            maxStrength = other.maxStrength;
            defaultStartColor = other.defaultStartColor;
            defaultEndColor = other.defaultEndColor;
        }
};

#endif