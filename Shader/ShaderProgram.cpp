#include "ShaderProgram.h"

ShaderProgram::ShaderProgram() {

}

ShaderProgram::ShaderProgram(unsigned int vertexId, unsigned int fragmentId) {
    this->shaderProgramId = glCreateProgram();
    glAttachShader(this->shaderProgramId, vertexId);
    glAttachShader(this->shaderProgramId, fragmentId);
    glLinkProgram(this->shaderProgramId);
}

void ShaderProgram::UseShaderProgram() {
    glUseProgram(this->shaderProgramId);
}

void ShaderProgram::DeleteLinkedShader(unsigned int vertexId, unsigned int fragmentId) {
    glDeleteShader(vertexId);
    glDeleteShader(fragmentId);
}

void ShaderProgram::SetInt(const std::string& uniformName, int val) {
    int uniformLocation = glGetUniformLocation(this->shaderProgramId, uniformName.c_str());
    glUniform1i(uniformLocation, val);
}

void ShaderProgram::SetMatrix4(const std::string& uniformName, glm::mat4 mat) {
    int uniformLocation = glGetUniformLocation(this->shaderProgramId, uniformName.c_str());
    glUniformMatrix4fv(uniformLocation, 1, GL_FALSE, glm::value_ptr(mat));
}