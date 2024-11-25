#ifndef _FLUID_SOLVER_H_
#define _FLUID_SOLVER_H_

#include <iostream>
#include <glad/glad.h> 
#include <GLFW/glfw3.h>
#include <vector>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "VectorField.h"

class FluidSolver {
    public:
        float timeStep;
        FluidSolver(float timeStep) {}

        VectorField Step(VectorField currentField);
};

#endif