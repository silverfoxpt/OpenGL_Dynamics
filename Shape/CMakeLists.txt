add_library(Shape
    TriangleDraw.cpp
    TriangleDraw.h
    TriangleColorDraw.cpp
    TriangleColorDraw.h
    TriangleColorTextureDraw.cpp
    TriangleColorTextureDraw.h
    SquareGridColorDraw.cpp
    SquareGridColorDraw.h
    LinesColorDraw.h
    LinesColorDraw.cpp
)

target_link_libraries(Shape PRIVATE glad glfw glm)