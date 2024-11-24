#ifndef _SHADER_PROGRAM_H_
#define _SHADER_PROGRAM_H_

#include <iostream>
#include <glad/glad.h> 
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

class ShaderProgram {
    public:
        unsigned int shaderProgramId;

        ShaderProgram();
        ShaderProgram(unsigned int vertexId, unsigned int fragmentId);

        void DeleteLinkedShader(unsigned int vertexId, unsigned int fragmentId);
        void UseShaderProgram();

        void SetInt(const std::string& uniformName, int val);
        void SetMatrix4(const std::string& uniformName, glm::mat4 mat);
};

#endif