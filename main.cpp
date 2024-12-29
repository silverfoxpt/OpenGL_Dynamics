#include <iostream>
#include <glad/glad.h> 
#include <GLFW/glfw3.h>
#include <cmath>
#include <thread>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "./Shader/VertexShader.h"
#include "./Shader/FragmentShader.h"
#include "./Shader/ShaderProgram.h"

#include "./Shape/TriangleDraw.h"
#include "./Shape/TriangleColorDraw.h"
#include "./Shape/TriangleColorTextureDraw.h"
#include "./Shape/SquareGridColorDraw.h"
#include "./Shape/LinesColorDraw.h"

#include "./Dynamics/VectorField.h"
#include "./Dynamics/ValueField.h"

// All texture implementation should be below here
#define STB_IMAGE_IMPLEMENTATION
#include "./Texture/Texture2D.h"

// Global variables for FPS calculation
auto lastTime = std::chrono::high_resolution_clock::now();
int frameCount = 0;

void PrintFPS() {
    // Current time
    auto currentTime = std::chrono::high_resolution_clock::now();
    
    // Elapsed time in seconds
    double elapsedTime = std::chrono::duration<double>(currentTime - lastTime).count();
    
    // Increment frame count
    frameCount++;
    
    // Print FPS every second
    if (elapsedTime >= 1.0) {
        std::cout << "FPS: " << frameCount / elapsedTime << std::endl;
        
        // Reset for the next interval
        frameCount = 0;
        lastTime = currentTime;
    }
}

GLFWwindow* window;
ShaderProgram colorOnlyShaderProgram;

// Values
int rows = 20, cols = 20;
float squareSize = 25;

VectorField vecField;
ValueField valField;

// Displaying
LinesColorDraw arrows;
SquareGridColorDraw colorSquareGrid;
std::vector<float> squareVertices; // For caching purposes

void Test() {
    stbi_set_flip_vertically_on_load(true); 
}

void ProcessDrawing() {
    // For Square Grid
    valField = ValueField(rows, cols, 0.0f);

    squareVertices  = SquareGridColorDraw::GenerateSampleGrid(rows, cols, squareSize);
    colorSquareGrid = SquareGridColorDraw(
        rows, 
        cols, 
        squareSize, 
        squareVertices, 
        valField.GenerateColorField(), 
        GL_DYNAMIC_DRAW
    );

    VertexShader colorOnlyVert = VertexShader(std::string(PROJECT_ROOT_DIR) + "/Shader/Vertex/colorOnly.vert");
    FragmentShader colorOnlyFrag = FragmentShader(std::string(PROJECT_ROOT_DIR) + "/Shader/Fragment/colorOnly.frag");

    colorOnlyShaderProgram = ShaderProgram(colorOnlyVert.shaderId, colorOnlyFrag.shaderId);
    colorOnlyShaderProgram.DeleteLinkedShader(colorOnlyVert.shaderId, colorOnlyFrag.shaderId);

    // For vector field & arrow
    vecField = VectorField(rows, cols, glm::vec2(10, -10));
    arrows = LinesColorDraw(
        vecField.GeneratePositionField(0, 0, squareSize), 
        vecField.GenerateColorField()
    );
}

void ProcessRendering() {
    // Get transformations
    glm::mat4 model = glm::mat4(1.0f); // world coords
    colorOnlyShaderProgram.SetMatrix4("model", model);
    
    glm::mat4 view = glm::mat4(1.0f); // camera matrix
    //view = glm::translate(view, glm::vec3(0.0f, 0.0f, -3.0f)); // Move camera backward 3 unit
    colorOnlyShaderProgram.SetMatrix4("view", view);

    glm::mat4 projection = glm::ortho(0.0f, 800.0f, -600.0f, 0.0f, 0.0f, 100.0f); //projection matrix - account for aspect ratio
    colorOnlyShaderProgram.SetMatrix4("projection", projection);

    // Draw square grid
    colorOnlyShaderProgram.UseShaderProgram();
    colorSquareGrid.UpdateColors(valField.GenerateColorField());
    colorSquareGrid.Draw(squareVertices.size() / 3);

    // Draw arrow grid
    arrows.UpdatePositions(vecField.GeneratePositionField(0, 0, squareSize));
    arrows.UpdateColors(vecField.GenerateColorField());
    arrows.Draw(2.0f);
}

/// @brief Resize window appropriately with glViewport, when user change window's size
/// @param window 
/// @param width 
/// @param height 
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
}  

/// @brief Let window check for inputs, and acts based on those inputs
/// @param window 
void processInput(GLFWwindow *window)
{
    // glfwGetKey - Check status of keys
    if(glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
}

int main() {
    // Initialization
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3); //version 4.3 then
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE); //using Core version

    // Create window with GLFW
    window = glfwCreateWindow(800, 600, "LearnOpenGL", NULL, NULL); //height, width, name
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();    
        return -1;
    }

    // Make the OpenGL context - very important!
    glfwMakeContextCurrent(window); 

    // Initialize GLAD - math library
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    } 

    // Create viewport - glViewport
    // First 2 set location of lower left corner of the rendered view
    // Last 2 set location of width and height of rendered view
    glViewport(0, 0, 800, 600);

    glDisable(GL_DEPTH_TEST);

    // Let user resize windows
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);  

    Test();
    ProcessDrawing();

    // Draw loop - Important
    while(!glfwWindowShouldClose(window))
    {
        // Input stuffs should be here
        processInput(window);

        // Drawing stuffs should be here
        // Clear color buffer
        glClearColor(0.2f, 0.3f, 0.3f, 1.0f); //state-setting!
        glClear(GL_COLOR_BUFFER_BIT); //state using!

        // Draw!
        ProcessRendering();
        
        PrintFPS();

        // Swap buffers
        // Many graphics system use two buffers. The inactive buffer serve as an intermediate 
        // for the machine to draw onto for the next frame, while the active buffer shows
        // the current frame. This way, there should be no flickers due to, e.g., race conditions.
        glfwSwapBuffers(window);
        
        // Check if any events are triggered - keyboards, mouse, etc...
        glfwPollEvents();    
    }

    // GLFW - Terminate when everything is done!
    glfwTerminate();
    return 0;
}
