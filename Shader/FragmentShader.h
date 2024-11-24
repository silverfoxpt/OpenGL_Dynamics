#ifndef _FRAGMENT_SHADER_H_
#define _FRAGMENT_SHADER_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <glad/glad.h> 
#include <GLFW/glfw3.h>

class FragmentShader {
    public:
        unsigned int shaderId;
        const char* shaderSource;

        FragmentShader(const std::string& filePath);
        void DeleteShader();

    private:
        int success;
        char infoLog[512];

        std::string ReadShaderFile(const std::string& filePath);
};

#endif