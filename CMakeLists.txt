cmake_minimum_required(VERSION 3.5.0)
project(OpenGLDynamics VERSION 0.1.0 LANGUAGES C CXX)

add_executable(OpenGLDynamics main.cpp)

add_definitions(-DPROJECT_ROOT_DIR=\"${PROJECT_SOURCE_DIR}\")

set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)

add_subdirectory("${PROJECT_SOURCE_DIR}/deps/glfw-3.4")
add_subdirectory("${PROJECT_SOURCE_DIR}/deps/glm")
add_subdirectory("${PROJECT_SOURCE_DIR}/Shader")
add_subdirectory("${PROJECT_SOURCE_DIR}/Shape")
add_subdirectory("${PROJECT_SOURCE_DIR}/Texture")
add_subdirectory("${PROJECT_SOURCE_DIR}/Dynamics")

# For OpenGL 4.6
add_subdirectory("${PROJECT_SOURCE_DIR}/deps/glad")
#include_directories("${PROJECT_SOURCE_DIR}/Image_Load")
#include_directories("${PROJECT_SOURCE_DIR}/Setup")

target_link_libraries(
    OpenGLDynamics
    glad
    glfw
    glm
    Shader
    Shape
    Texture
    Dynamics
)

