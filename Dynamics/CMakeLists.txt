add_library(Dynamics
    VectorField.h
    VectorField.cpp
)

target_link_libraries(Dynamics PRIVATE glad glfw glm)