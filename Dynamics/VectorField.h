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
        float maxStrength = 10.0;
        glm::vec3 defaultStartColor = glm::vec3(1.0f, 0.0f, 0.0f); //Red
        glm::vec3 defaultEndColor = glm::vec3(1.0f, 0.0f, 0.0f); //green

        VectorField() { }
        VectorField(int rows, int cols, glm::vec2 initialVec);

        glm::vec2 GetVector(int rowIdx, int colIdx);
        glm::vec2 GetUnitVector(int rowIdx, int colIdx);
        float GetStrengthVector(int rowIdx, int colIdx);

        std::vector<float> GeneratePositionField(float startX, float startY, float length);
        std::vector<float> GenerateColorField();
};

#endif