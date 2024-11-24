add_library(Texture
    Texture2D.cpp
    Texture2D.h
)

target_link_libraries(Texture PRIVATE glad glfw glm)