#ifndef _FLUID_SOLVER_ARAKAWAC_H_
#define _FLUID_SOLVER_ARAKAWAC_H_

#include <iostream>
#include <glad/glad.h> 
#include <GLFW/glfw3.h>
#include <vector>
#include <thread>
#include <cmath>
#include <random>
#include <omp.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "VectorField.h"
#include "ValueField.h"

class FluidSolverC {
    public:
        struct BlobInflow {
            int lifeLeft;
            int radius;
            float strength;
            glm::vec2 direction;

            int jCenter;                             // horizontal center
            int currentLayer;                        // current row layer being injected
            std::vector<std::vector<float>> layers;  // precomputed circular blob slices
        };

        float timeStep;
        float initialDensity = 0.0f;
        float diffusionRate = 0.01f;
        glm::vec2 initialVelocity = glm::vec2(0.0f, 0.0f);

        int rows, cols;
        int circleRadius = 3; // radius of solid circle obstacles
        int startingRow = 160; // row to start injecting blobs
        float currentAngle = 0;
        int gaussSeidelIterations = 20;

        ValueField currentDensity, nextDensity;

        // Velocity fields (staggered)
        ValueField uVel, nextUVel; // x-direction velocity: size (rows, cols+1)
        ValueField vVel, nextVVel; // y-direction velocity: size (rows+1, cols)

        ValueField divergence;
        ValueField pressure;
        ValueField obstacleMask;

        FluidSolverC(int rows, int cols);

        void Diffusion();
        void Advection();
        void ClearDivergence();

        void SetReflectiveBoundary();

        void SetVonNeumannBoundary();
        void SetSpongeBoundary(bool disperseVelocityOnly);
        void ApplyRadiativeBoundary(ValueField& field, float waveSpeed);

        void InjectRandomEddiesAtTop(float chance);
        void InjectInstantBlobAtTop(float strength, int width, float density);

        std::vector<BlobInflow> activeBlobs;
        void InjectBlob(BlobInflow& blob);
        void MaybeSpawnBlob();

        void AddSolidCircle(int centerI, int centerJ, int radius);

        void Step();

        // Interpolate velocity at cell centers
        glm::vec2 GetVelocityAtCenter(int i, int j) const;

        // Build vertex positions for the velocity arrows
        std::vector<float> GenerateArrowPositions(float startX,
            float startY,
            float cellSize,
            float maxStrength = 60.0f) const;

        // Build per-vertex colours (red → green ramp)
        std::vector<float> GenerateArrowColours(float maxStrength = 60.0f,
        const glm::vec3& slowColour = {1,0,0},
        const glm::vec3& fastColour = {0,1,0}) const;
};

#endif