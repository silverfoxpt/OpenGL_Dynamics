add_library(Shader
    VertexShader.cpp
    VertexShader.h
    FragmentShader.cpp
    FragmentShader.h
    ShaderProgram.cpp
    ShaderProgram.h
)

target_link_libraries(Shader PRIVATE glad glfw glm)