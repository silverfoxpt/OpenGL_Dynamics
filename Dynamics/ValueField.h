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
        float maxValue = 1.0;
        glm::vec3 defaultMaxColor = glm::vec3(1.0f, 1.0f, 1.0f); //White

        ValueField() { }
        ValueField(int rows, int cols, float initialVal);

        float GetValue(int rowIdx, int colIdx);
        void SetValue(int rowIdx, int colIdx, float newVal);
        
        std::vector<float> GenerateColorField();
};

#endif