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
#include "ValueField.h"

class FluidSolver {
    public:
        float timeStep;
        float initialDensity = 0.0f;
        float diffusionRate = 0.05f;
        glm::vec2 initialVelocity = glm::vec2(0.0f, 0.0f);

        int rows, cols;

        ValueField currentDensity, nextDensity;
        VectorField currentVelocity, nextVelocity;

        FluidSolver(int rows, int cols);

        void Diffusion();
        void Advection();
        void ClearDivergence();
        void AddVorticityConfinement();

        void SetReflectiveBoundary();
        void Step();
};

#endif