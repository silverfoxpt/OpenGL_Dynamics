#include "FluidSolver.h"

FluidSolver::FluidSolver(int rows, int cols) {
    // Initialize
    this->rows = rows; this->cols = cols; this->timeStep = timeStep;
    currentDensity = nextDensity = ValueField(rows, cols, initialDensity);
    currentVelocity = nextVelocity = VectorField(rows, cols, initialVelocity);
}

void FluidSolver::SetReflectiveBoundary() {
    // Reflective boundary conditions for the velocity field (no penetration)
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            if (i == 0) { // Top row
                auto velocity = currentVelocity.GetVector(i+1, j);
                currentVelocity.SetVector(i, j, glm::vec2(velocity.x, -velocity.y)); // Reverse y
                currentDensity.SetValue(i, j, currentDensity.GetValue(i + 1, j));   // Mirror density
            }
            if (i == rows - 1) { // Bottom row
                auto velocity = currentVelocity.GetVector(i-1, j);
                currentVelocity.SetVector(i, j, glm::vec2(velocity.x, -velocity.y)); // Reverse y
                currentDensity.SetValue(i, j, currentDensity.GetValue(i - 1, j));   // Mirror density
            }
            if (j == 0) { // Left column
                auto velocity = currentVelocity.GetVector(i, j+1);
                currentVelocity.SetVector(i, j, glm::vec2(-velocity.x, velocity.y)); // Reverse x
                currentDensity.SetValue(i, j, currentDensity.GetValue(i, j + 1));   // Mirror density
            }
            if (j == cols - 1) { // Right column
                auto velocity = currentVelocity.GetVector(i, j-1);
                currentVelocity.SetVector(i, j, glm::vec2(-velocity.x, velocity.y)); // Reverse x
                currentDensity.SetValue(i, j, currentDensity.GetValue(i, j - 1));   // Mirror density
            }
        }
    }
}

void FluidSolver::Diffusion() {
    SetReflectiveBoundary();

    // Gauss-Seidel parameters
    int maxIterations = 3; // Number of iterations for convergence
    float alpha = diffusionRate * timeStep;

    // Function to process a 2x2 square region
    auto processSquare = [&](int startRow, int startCol, int endRow, int endCol) {
        for (int i = startRow; i < endRow; ++i) {
            for (int j = startCol; j < endCol; ++j) {
                if (i >= 1 && i < rows - 1 && j >= 1 && j < cols - 1) {
                    // Accumulate values of valid neighbors
                    float sum = 0.0f;
                    int count = 0;

                    // Check neighbors
                    if (i > 0) { sum += currentDensity.GetValue(i - 1, j); ++count; }
                    if (i < rows - 1) { sum += currentDensity.GetValue(i + 1, j); ++count; }
                    if (j > 0) { sum += currentDensity.GetValue(i, j - 1); ++count; }
                    if (j < cols - 1) { sum += currentDensity.GetValue(i, j + 1); ++count; }

                    // Average neighbors
                    float avg = sum / count;

                    // Update cell value
                    float newVal = (currentDensity.GetValue(i, j) + alpha * avg) / (1 + alpha);
                    nextDensity.SetValue(i, j, newVal);
                }
            }
        }
    };

    for (int iter = 0; iter < maxIterations; ++iter) {
        const int numThreads = 4; // Fixed number of threads
        std::vector<std::thread> threads;

        // Assign 2x2 blocks to threads
        int blockHeight = (rows - 2) / 2; // Divide grid into two vertical parts
        int blockWidth = (cols - 2) / 2;  // Divide grid into two horizontal parts

        threads.emplace_back(processSquare, 1, 1, 1 + blockHeight, 1 + blockWidth);                     // Top-left
        threads.emplace_back(processSquare, 1, 1 + blockWidth, 1 + blockHeight, cols - 1);             // Top-right
        threads.emplace_back(processSquare, 1 + blockHeight, 1, rows - 1, 1 + blockWidth);             // Bottom-left
        threads.emplace_back(processSquare, 1 + blockHeight, 1 + blockWidth, rows - 1, cols - 1);      // Bottom-right

        // Join all threads
        for (auto& thread : threads) {
            thread.join();
        }
    }

    currentDensity.SwapWith(nextDensity);
}

void FluidSolver::Advection() {
    SetReflectiveBoundary();

    // Parameters for advection
    float dt = timeStep;
    float gridSpacing = 1.0f; // Assuming uniform grid spacing of 1 unit.

    // Function to process a 2x2 block
    auto processBlock = [&](int startRow, int startCol, int endRow, int endCol) {
        for (int i = startRow; i < endRow; ++i) {
            for (int j = startCol; j < endCol; ++j) {
                if (i >= 1 && i < rows - 1 && j >= 1 && j < cols - 1) {
                    // Compute the backward position of the current cell
                    float x = j - currentVelocity.GetVector(i, j).x * dt / gridSpacing;
                    float y = i - currentVelocity.GetVector(i, j).y * dt / gridSpacing;

                    // Use bilinear interpolation to get density at the backward position
                    float dx = x - floor(x);  // Fractional part of x
                    float dy = y - floor(y);  // Fractional part of y

                    int x0 = static_cast<int>(floor(x));
                    int y0 = static_cast<int>(floor(y));
                    int x1 = x0 + 1;
                    int y1 = y0 + 1;

                    // Ensure coordinates are within bounds
                    x0 = std::max(0, std::min(x0, cols - 1));
                    x1 = std::max(0, std::min(x1, cols - 1));
                    y0 = std::max(0, std::min(y0, rows - 1));
                    y1 = std::max(0, std::min(y1, rows - 1));

                    // Perform bilinear interpolation for density
                    float topLeft = currentDensity.GetValue(y0, x0);
                    float topRight = currentDensity.GetValue(y0, x1);
                    float bottomLeft = currentDensity.GetValue(y1, x0);
                    float bottomRight = currentDensity.GetValue(y1, x1);

                    float interpolatedDensity =
                        (1 - dx) * (1 - dy) * topLeft +
                        dx * (1 - dy) * topRight +
                        (1 - dx) * dy * bottomLeft +
                        dx * dy * bottomRight;

                    // Set the interpolated density value at the new location
                    nextDensity.SetValue(i, j, interpolatedDensity);

                    // Perform bilinear interpolation for velocity components
                    glm::vec2 topLeftVel = currentVelocity.GetVector(y0, x0);
                    glm::vec2 topRightVel = currentVelocity.GetVector(y0, x1);
                    glm::vec2 bottomLeftVel = currentVelocity.GetVector(y1, x0);
                    glm::vec2 bottomRightVel = currentVelocity.GetVector(y1, x1);

                    float interpolatedVelX =
                        (1 - dx) * (1 - dy) * topLeftVel.x +
                        dx * (1 - dy) * topRightVel.x +
                        (1 - dx) * dy * bottomLeftVel.x +
                        dx * dy * bottomRightVel.x;

                    float interpolatedVelY =
                        (1 - dx) * (1 - dy) * topLeftVel.y +
                        dx * (1 - dy) * topRightVel.y +
                        (1 - dx) * dy * bottomLeftVel.y +
                        dx * dy * bottomRightVel.y;

                    // Set the interpolated velocity value at the new location
                    nextVelocity.SetVector(i, j, {interpolatedVelX, interpolatedVelY});
                }
            }
        }
    };

    const int numThreads = 4; // Fixed number of threads
    std::vector<std::thread> threads;

    // Assign 2x2 blocks to threads
    int blockHeight = (rows - 2) / 2; // Divide grid into two vertical parts
    int blockWidth = (cols - 2) / 2;  // Divide grid into two horizontal parts

    threads.emplace_back(processBlock, 1, 1, 1 + blockHeight, 1 + blockWidth);                     // Top-left
    threads.emplace_back(processBlock, 1, 1 + blockWidth, 1 + blockHeight, cols - 1);             // Top-right
    threads.emplace_back(processBlock, 1 + blockHeight, 1, rows - 1, 1 + blockWidth);             // Bottom-left
    threads.emplace_back(processBlock, 1 + blockHeight, 1 + blockWidth, rows - 1, cols - 1);      // Bottom-right

    // Join all threads
    for (auto& thread : threads) {
        thread.join();
    }

    // Update the current density and velocity to the next ones after advection
    currentDensity.SwapWith(nextDensity);
    currentVelocity.SwapWith(nextVelocity);
}

void FluidSolver::ClearDivergence() {
    SetReflectiveBoundary();

    // Parameters for the pressure solver
    const int maxIterations = 3;  // Iterations for pressure solving
    const float overRelaxation = 1.9f;  // Over-relaxation parameter for faster convergence
    ValueField pressure(rows, cols, 0.0f);  // Pressure field to solve for

    // First calculate divergence at each cell using central differences
    ValueField divergence(rows, cols, 0.0f);
    
    // Function to process a 2x2 block for divergence calculation
    auto processDivergenceBlock = [&](int startRow, int startCol, int endRow, int endCol) {
        for (int i = startRow; i < endRow; ++i) {
            for (int j = startCol; j < endCol; ++j) {
                if (i >= 1 && i < rows - 1 && j >= 1 && j < cols - 1) {
                    // Calculate divergence using the formula from the image
                    float vxRight = currentVelocity.GetVector(i, j + 1).x;
                    float vxLeft = currentVelocity.GetVector(i, j - 1).x;
                    float vyUp = currentVelocity.GetVector(i - 1, j).y;
                    float vyDown = currentVelocity.GetVector(i + 1, j).y;

                    divergence.SetValue(i, j, (vxRight - vxLeft + vyDown - vyUp) / 2.0f);
                }
            }
        }
    };

    // Divide the grid into 4 blocks
    int blockHeight = (rows - 2) / 2;
    int blockWidth = (cols - 2) / 2;
    std::vector<std::thread> threads;

    threads.emplace_back(processDivergenceBlock, 1, 1, 1 + blockHeight, 1 + blockWidth);        // Top-left block
    threads.emplace_back(processDivergenceBlock, 1, 1 + blockWidth, 1 + blockHeight, cols - 1);    // Top-right block
    threads.emplace_back(processDivergenceBlock, 1 + blockHeight, 1, rows - 1, 1 + blockWidth);    // Bottom-left block
    threads.emplace_back(processDivergenceBlock, 1 + blockHeight, 1 + blockWidth, rows - 1, cols - 1); // Bottom-right block

    // Join all threads
    for (auto& thread : threads) {
        thread.join();
    }

    // Solve the pressure Poisson equation using Gauss-Seidel relaxation
    // Function to process a 2x2 block for pressure solving
    auto processPressureBlock = [&](int startRow, int startCol, int endRow, int endCol) {
        for (int iter = 0; iter < maxIterations; ++iter) {
            for (int i = startRow; i < endRow; ++i) {
                for (int j = startCol; j < endCol; ++j) {
                    if (i >= 1 && i < rows - 1 && j >= 1 && j < cols - 1) {
                        float pLeft = pressure.GetValue(i, j - 1);
                        float pRight = pressure.GetValue(i, j + 1);
                        float pUp = pressure.GetValue(i - 1, j);
                        float pDown = pressure.GetValue(i + 1, j);

                        // Calculate new pressure using over-relaxation
                        float pCenter = pressure.GetValue(i, j);
                        float newPressure = (pLeft + pRight + pUp + pDown - divergence.GetValue(i, j)) / 4.0f;
                        pressure.SetValue(i, j, pCenter + overRelaxation * (newPressure - pCenter));
                    }
                }
            }
        }
    };

    // Divide the grid into 4 blocks for pressure solving
    threads.clear();
    threads.emplace_back(processPressureBlock, 1, 1, 1 + blockHeight, 1 + blockWidth);           // Top-left block
    threads.emplace_back(processPressureBlock, 1, 1 + blockWidth, 1 + blockHeight, cols - 1);     // Top-right block
    threads.emplace_back(processPressureBlock, 1 + blockHeight, 1, rows - 1, 1 + blockWidth);     // Bottom-left block
    threads.emplace_back(processPressureBlock, 1 + blockHeight, 1 + blockWidth, rows - 1, cols - 1); // Bottom-right block

    // Join all threads
    for (auto& thread : threads) {
        thread.join();
    }

    // Subtract pressure gradient from velocity field to make it divergence-free
    // Function to process a 2x2 block for velocity update
    auto processVelocityBlock = [&](int startRow, int startCol, int endRow, int endCol) {
        for (int i = startRow; i < endRow; ++i) {
            for (int j = startCol; j < endCol; ++j) {
                if (i >= 1 && i < rows - 1 && j >= 1 && j < cols - 1) {
                    // Calculate pressure gradient using central differences
                    float gradPx = (pressure.GetValue(i, j + 1) - pressure.GetValue(i, j - 1)) / 2.0f;
                    float gradPy = (pressure.GetValue(i + 1, j) - pressure.GetValue(i - 1, j)) / 2.0f;

                    // Subtract gradient from velocity field
                    glm::vec2 currentVel = currentVelocity.GetVector(i, j);
                    currentVelocity.SetVector(i, j, {
                        currentVel.x - gradPx,
                        currentVel.y - gradPy
                    });
                }
            }
        }
    };

    // Divide the grid into 4 blocks for velocity update
    threads.clear();
    threads.emplace_back(processVelocityBlock, 1, 1, 1 + blockHeight, 1 + blockWidth);           // Top-left block
    threads.emplace_back(processVelocityBlock, 1, 1 + blockWidth, 1 + blockHeight, cols - 1);     // Top-right block
    threads.emplace_back(processVelocityBlock, 1 + blockHeight, 1, rows - 1, 1 + blockWidth);     // Bottom-left block
    threads.emplace_back(processVelocityBlock, 1 + blockHeight, 1 + blockWidth, rows - 1, cols - 1); // Bottom-right block

    // Join all threads
    for (auto& thread : threads) {
        thread.join();
    }

    //SetReflectiveBoundary();
}

void FluidSolver::Step() {
    // Add upward velocity to create movement
    float baseVelocity = currentVelocity.maxStrength; // Adjust this value to control speed
    glm::vec2 sourceVelocity(baseVelocity, 0.0f); // Negative y means upward

    int offset = 40;
    for (int i = rows / 2 - offset; i < rows / 2 + offset; i++) {
        for (int j = cols / 2 - offset; j < cols / 2 + offset; j++) {
            this->currentDensity.SetValue(i, j, 1.0f);
        }
    }

    for (int j = cols / 2 - offset; j < cols - 1; j++) {
        for (int i = rows / 2 - offset; i < rows / 2 + offset; i++) {
            this->currentVelocity.SetVector(rows / 2, j, sourceVelocity);
        }
    }

    this->Diffusion();
    this->Advection();
    this->ClearDivergence();
}

