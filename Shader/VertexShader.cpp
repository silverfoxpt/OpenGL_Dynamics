#include "VertexShader.h"

VertexShader::VertexShader(const std::string& filePath) {
    std::string fileContent = this->ReadShaderFile(filePath);
    this->shaderSource = fileContent.c_str();

    this->shaderId = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(this->shaderId, 1, &this->shaderSource, NULL);
    glCompileShader(this->shaderId);

    // Check for error
    glGetShaderiv(this->shaderId, GL_COMPILE_STATUS, &this->success);
    if (!success) {
        glGetShaderInfoLog(this->shaderId, 512, NULL, this->infoLog);
        throw new std::runtime_error("Shader compile failed: " + std::string(this->infoLog));
    }
}

void VertexShader::DeleteShader() {
    glDeleteShader(this->shaderId);
}

std::string VertexShader::ReadShaderFile(const std::string& filePath) {
    std::ifstream shaderFile(filePath);   // Open the file
    if (!shaderFile.is_open()) {
        throw std::invalid_argument("File does not exist");
    }

    std::stringstream shaderStream;
    shaderStream << shaderFile.rdbuf();   // Read file into string stream
    shaderFile.close();                   // Close file

    return shaderStream.str();            // Return the shader code as a string
}