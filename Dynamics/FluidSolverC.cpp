#include "FluidSolverC.h"

FluidSolverC::FluidSolverC(int rows, int cols) {
    this->rows = rows;
    this->cols = cols;
    this->timeStep = 0.1f;  // default value or configurable externally

    // Initialize scalar fields on the cell centers
    currentDensity = nextDensity = ValueField(rows, cols, initialDensity);
    pressure = divergence = ValueField(rows, cols, 0.0f);

    // Initialize velocity components staggered in C-grid style:
    // u (horizontal) lives on vertical edges => (rows, cols+1)
    // v (vertical) lives on horizontal edges => (rows+1, cols)
    uVel = nextUVel = ValueField(rows, cols + 1, initialVelocity.x); // u-component
    vVel = nextVVel = ValueField(rows + 1, cols, initialVelocity.y); // v-component
}

glm::vec2 FluidSolverC::GetVelocityAtCenter(int i, int j) const
{
    float u_c = 0.5f * (uVel.valueField[i][j] + uVel.valueField[i][j + 1]);
    float v_c = 0.5f * (vVel.valueField[i][j] + vVel.valueField[i + 1][j]);
    return glm::vec2(u_c, v_c);
}

// Newmann OBC
void FluidSolverC::SetVonNeumannBoundary()
{
    // Density: copy adjacent values (zero normal gradient)
    for (int j = 0; j < cols; ++j) {
        currentDensity.SetValue(0,        j, currentDensity.GetValue(1,        j)); // top
        currentDensity.SetValue(rows - 1, j, currentDensity.GetValue(rows - 2, j)); // bottom
    }
    for (int i = 0; i < rows; ++i) {
        currentDensity.SetValue(i, 0,        currentDensity.GetValue(i, 1));        // left
        currentDensity.SetValue(i, cols - 1, currentDensity.GetValue(i, cols - 2)); // right
    }

    // u velocity: extrapolate or zero-gradient (free outflow)
    for (int i = 0; i < rows; ++i) {
        uVel.SetValue(i, 0,     uVel.GetValue(i, 1));          // left
        uVel.SetValue(i, cols,  uVel.GetValue(i, cols - 1));   // right
    }

    for (int j = 0; j <= cols; ++j) {
        uVel.SetValue(0,        j, uVel.GetValue(1,        j)); // top
        uVel.SetValue(rows - 1, j, uVel.GetValue(rows - 2, j)); // bottom
    }

    // v velocity: extrapolate or zero-gradient
    for (int j = 0; j < cols; ++j) {
        vVel.SetValue(0,     j, vVel.GetValue(1,     j));     // top
        vVel.SetValue(rows,  j, vVel.GetValue(rows - 1, j));  // bottom
    }

    for (int i = 0; i <= rows; ++i) {
        vVel.SetValue(i, 0,        vVel.GetValue(i, 1));        // left
        vVel.SetValue(i, cols - 1, vVel.GetValue(i, cols - 2)); // right
    }
}

// Sponge OBC
void FluidSolverC::SetSpongeBoundary(bool disperseVelocityOnly = true)
{
    const int   spongeThickness     = 20;
    const float spongeFadeExponent  = 0.09f;

    auto fade = [&](float alpha)
    {
        //return std::max(0.99f, std::pow(1.0f - alpha, spongeFadeExponent));
        return std::pow(1.0f - alpha, spongeFadeExponent);
    };

    // === Density: zero-gradient, optionally with sponge damping ===
    for (int j = 0; j < cols; ++j) {
        for (int i = 0; i < spongeThickness; ++i) {
            float alpha = float(i) / spongeThickness;
            float damp  = fade(alpha);
            float valT  = currentDensity.GetValue(i + 1, j);
            float valB  = currentDensity.GetValue(rows - 2 - i, j);

            if (disperseVelocityOnly) {
                currentDensity.SetValue(i,            j, valT);  // zero-gradient only
                currentDensity.SetValue(rows - 1 - i, j, valB);
            } else {
                currentDensity.SetValue(i,            j, valT * damp);  // sponge
                currentDensity.SetValue(rows - 1 - i, j, valB * damp);
            }
        }
    }

    for (int i = 0; i < rows; ++i) {
        currentDensity.SetValue(i, 0,        currentDensity.GetValue(i, 1));
        currentDensity.SetValue(i, cols - 1, currentDensity.GetValue(i, cols - 2));
    }

    // === u velocity: zero-gradient with sponge ===
    for (int j = 0; j <= cols; ++j) {
        for (int i = 0; i < spongeThickness; ++i) {
            float alpha = float(i) / spongeThickness;
            float damp  = fade(alpha);
            float valT  = uVel.GetValue(i + 1, j);
            float valB  = uVel.GetValue(rows - 2 - i, j);

            uVel.SetValue(i,             j, valT * damp);
            uVel.SetValue(rows - 1 - i,  j, valB * damp);
        }
    }

    for (int i = 0; i < rows; ++i) {
        uVel.SetValue(i, 0,     uVel.GetValue(i, 1));
        uVel.SetValue(i, cols,  uVel.GetValue(i, cols - 1));
    }

    // === v velocity: zero-gradient with sponge ===
    for (int j = 0; j < cols; ++j) {
        for (int i = 0; i < spongeThickness; ++i) {
            float alpha = float(i) / spongeThickness;
            float damp  = fade(alpha);

            float valT = vVel.GetValue(i + 1, j);           // top
            float valB = vVel.GetValue(rows - 2 - i, j);    // bottom

            // Top: only damp upward (positive) velocities
            if (valT > 0.0f)
                vVel.SetValue(i, j, valT * damp);
            else
                vVel.SetValue(i, j, valT);

            // Bottom: only damp downward (negative) velocities
            if (valB < 0.0f)
                vVel.SetValue(rows - i, j, valB * damp);
            else
                vVel.SetValue(rows - i, j, valB);
        }
    }
}

void FluidSolverC::SetReflectiveBoundary()
{
    for (int j = 0; j < cols; ++j) {                 // top / bottom
        currentDensity.SetValue(0,        j, currentDensity.GetValue(1,        j));
        currentDensity.SetValue(rows - 1, j, currentDensity.GetValue(rows - 2, j));
    }
    for (int i = 0; i < rows; ++i) {                 // left / right
        currentDensity.SetValue(i, 0,        currentDensity.GetValue(i, 1));
        currentDensity.SetValue(i, cols - 1, currentDensity.GetValue(i, cols - 2));
    }

    for (int i = 0; i < rows; ++i) {       // left & right walls (normal = x)
        uVel.SetValue(i, 0,     0.0f);     // left wall  (face index 0)
        uVel.SetValue(i, cols,  0.0f);     // right wall (face index cols)
    }
    for (int j = 0; j <= cols; ++j) {      // top & bottom (tangential)
        uVel.SetValue(0,        j, uVel.GetValue(1,        j));          // top
        uVel.SetValue(rows - 1, j, uVel.GetValue(rows - 2, j));          // bottom
    }

    for (int j = 0; j < cols; ++j) {       // top & bottom walls (normal = y)
        vVel.SetValue(0,     j, 0.0f);     // top wall    (face index 0)
        vVel.SetValue(rows,  j, 0.0f);     // bottom wall (face index rows)
    }
    for (int i = 0; i <= rows; ++i) {      // left & right (tangential)
        vVel.SetValue(i, 0,        vVel.GetValue(i, 1));                  // left
        vVel.SetValue(i, cols - 1, vVel.GetValue(i, cols - 2));           // right
    }
}

void FluidSolverC::ApplyRadiativeBoundary(ValueField& field, float waveSpeed)
{
    const float dt = timeStep;
    const float dx = 1.0f;  // assuming uniform grid spacing

    int numCols = field.cols;

    // Apply to the top boundary (row = 0)
    for (int j = 0; j < numCols; ++j)
    {
        // Use two interior values: row 1 and 2
        float phi0 = field.GetValue(1, j);  // just inside domain
        float phi1 = field.GetValue(2, j);  // one more inside

        // Radiative condition (Sommerfeld-type)
        float phiOut = phi0 - waveSpeed * dt * (phi0 - phi1) / dx;

        // Set top boundary (ghost cell)
        field.SetValue(0, j, phiOut);
    }
}

void FluidSolverC::Diffusion()
{
    //SetReflectiveBoundary();
    //SetOpenBoundary();

    const float alpha = diffusionRate * timeStep;
    const int blockH = (rows - 2) / 2;          // interior rows per stripe
    const int blockW = (cols - 2) / 4;          // interior cols per stripe

    auto relaxBlock = [&](int r0,int c0,int r1,int c1)
    {
        #pragma omp parallel for collapse(2)
        for (int i = r0; i < r1; ++i)
            for (int j = c0; j < c1; ++j)
            {
                if (i < 1 || i >= rows-1 || j < 1 || j >= cols-1) continue;

                float sum =
                      currentDensity.GetValue(i - 1, j)   // north
                    + currentDensity.GetValue(i + 1, j)   // south
                    + currentDensity.GetValue(i, j - 1)   // west
                    + currentDensity.GetValue(i, j + 1);  // east

                float avg    = 0.25f * sum;
                float newVal = (currentDensity.GetValue(i, j) + alpha * avg)
                               / (1.0f + alpha);

                nextDensity.SetValue(i, j, newVal);
            }
    };

    for (int sweep = 0; sweep < gaussSeidelIterations; ++sweep)
    {
        std::vector<std::thread> th;
        for (int bi = 0; bi < 2; ++bi)
            for (int bj = 0; bj < 4; ++bj)
            {
                int r0 = 1 +  bi      * blockH;
                int r1 = 1 + (bi + 1) * blockH;
                int c0 = 1 +  bj      * blockW;
                int c1 = 1 + (bj + 1) * blockW;
                th.emplace_back(relaxBlock, r0, c0, r1, c1);
            }
        for (auto& t : th) t.join();
    }

    currentDensity.SwapWith(nextDensity);
}

void FluidSolverC::Advection()
{
    //SetReflectiveBoundary();
    //SetVonNeumannBoundary();
    //SetSpongeBoundary();
    ApplyRadiativeBoundary(vVel, 25.0f);  // Top edge of vertical velocity
    //ApplyRadiativeBoundary(currentDensity, 0.5f);  // Dye propagation

    const float dt    = timeStep;
    const float h     = 1.0f;
    const float inv_h = 1.0f / h;

    float maxVel = 0.0f;
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            maxVel = std::max(maxVel,
                              glm::length(GetVelocityAtCenter(i, j)));

    if (maxVel * dt * inv_h > 1.0f)
        std::cerr << "[CFL WARNING]  max(|v|)*dt = "
                  << maxVel * dt * inv_h << " > 1\n";

    auto clampInt = [](int v, int lo, int hi)
    { return std::max(lo, std::min(v, hi)); };

    const int blockH = (rows - 2) / 2;        // interior rows per stripe
    const int blockW = (cols - 2) / 4;        // interior cols per stripe

    auto advectDensityBlock = [&](int r0,int c0,int r1,int c1)
    {
        for (int i = r0; i < r1; ++i)
            for (int j = c0; j < c1; ++j)
            {
                if (i < 1 || i >= rows-1 || j < 1 || j >= cols-1) continue;

                glm::vec2 vel = GetVelocityAtCenter(i, j);
                float x = j - vel.x * dt * inv_h;
                float y = i - vel.y * dt * inv_h;

                int   x0 = int(std::floor(x)),  y0 = int(std::floor(y));
                float dx = x - x0,              dy = y - y0;

                x0 = clampInt(x0, 0, cols - 1);
                int x1 = clampInt(x0 + 1, 0, cols - 1);
                y0 = clampInt(y0, 0, rows - 1);
                int y1 = clampInt(y0 + 1, 0, rows - 1);

                float ρ =
                    (1 - dx) * (1 - dy) * currentDensity.GetValue(y0, x0) +
                    dx       * (1 - dy) * currentDensity.GetValue(y0, x1) +
                    (1 - dx) * dy       * currentDensity.GetValue(y1, x0) +
                    dx       * dy       * currentDensity.GetValue(y1, x1);

                nextDensity.SetValue(i, j, ρ);
            }
    };

    {
        std::vector<std::thread> th;
        for (int bi = 0; bi < 2; ++bi)
            for (int bj = 0; bj < 4; ++bj)
            {
                int r0 = 1 +  bi      * blockH;
                int r1 = 1 + (bi + 1) * blockH;
                int c0 = 1 +  bj      * blockW;
                int c1 = 1 + (bj + 1) * blockW;
                th.emplace_back(advectDensityBlock, r0, c0, r1, c1);
            }
        for (auto& t : th) t.join();
    }

    auto advectUBlock = [&](int r0,int c0,int r1,int c1)
    {
        for (int i = r0; i < r1; ++i)
            for (int j = c0; j < c1; ++j)
            {
                if (j < 1 || j >= cols) continue;           // walls

                float u = uVel.GetValue(i, j);
                float v = 0.5f * ( vVel.GetValue(i,     j-1)
                                 + vVel.GetValue(i + 1, j-1) );

                glm::vec2 vel(u, v);

                float x = j - vel.x * dt * inv_h;
                float y = i - vel.y * dt * inv_h;

                int   x0 = int(std::floor(x)),  y0 = int(std::floor(y));
                float dx = x - x0,              dy = y - y0;

                x0 = clampInt(x0, 0, cols);
                int x1 = clampInt(x0 + 1, 0, cols);
                y0 = clampInt(y0, 0, rows - 1);
                int y1 = clampInt(y0 + 1, 0, rows - 1);

                float uNew =
                    (1 - dx) * (1 - dy) * uVel.GetValue(y0, x0) +
                    dx       * (1 - dy) * uVel.GetValue(y0, x1) +
                    (1 - dx) * dy       * uVel.GetValue(y1, x0) +
                    dx       * dy       * uVel.GetValue(y1, x1);

                nextUVel.SetValue(i, j, uNew);
            }
    };

    {
        std::vector<std::thread> th;
        for (int bi = 0; bi < 2; ++bi)
            for (int bj = 0; bj < 4; ++bj)
            {
                int r0 =      bi      * blockH;                   // include halo rows
                int r1 = (bi == 1) ? rows : r0 + blockH;
                int c0 =  bj      * blockW;
                int c1 = (bj == 3) ? cols+1 : c0 + blockW;
                th.emplace_back(advectUBlock, r0, c0, r1, c1);
            }
        for (auto& t : th) t.join();
    }

    auto advectVBlock = [&](int r0,int c0,int r1,int c1)
    {
        for (int i = r0; i < r1; ++i)
        {
            if (i < 1 || i >= rows) continue;                 // walls
            for (int j = c0; j < c1; ++j)
            {
                float v = vVel.GetValue(i, j);
                float u = 0.5f * ( uVel.GetValue(i-1, j + 1)
                                 + uVel.GetValue(i-1, j) );

                glm::vec2 vel(u, v);

                float x = j - vel.x * dt * inv_h;
                float y = i - vel.y * dt * inv_h;

                int   x0 = int(std::floor(x)),  y0 = int(std::floor(y));
                float dx = x - x0,              dy = y - y0;

                x0 = clampInt(x0, 0, cols - 1);
                int x1 = clampInt(x0 + 1, 0, cols - 1);
                y0 = clampInt(y0, 0, rows);
                int y1 = clampInt(y0 + 1, 0, rows);

                float vNew =
                    (1 - dx) * (1 - dy) * vVel.GetValue(y0, x0) +
                    dx       * (1 - dy) * vVel.GetValue(y0, x1) +
                    (1 - dx) * dy       * vVel.GetValue(y1, x0) +
                    dx       * dy       * vVel.GetValue(y1, x1);

                nextVVel.SetValue(i, j, vNew);
            }
        }
    };

    {
        std::vector<std::thread> th;
        for (int bi = 0; bi < 2; ++bi)
            for (int bj = 0; bj < 4; ++bj)
            {
                int r0 =  bi      * blockH;
                int r1 = (bi == 1) ? rows+1 : r0 + blockH;
                int c0 = 1 +  bj      * blockW;      // same interior span as density
                int c1 = 1 + (bj + 1) * blockW;
                th.emplace_back(advectVBlock, r0, c0, r1, c1);
            }
        for (auto& t : th) t.join();
    }

    currentDensity.SwapWith(nextDensity);
    uVel.SwapWith(nextUVel);
    vVel.SwapWith(nextVVel);

    //SetReflectiveBoundary();
    //SetOpenBoundary();
}

void FluidSolverC::ClearDivergence()
{
    const float dx      = 1.0f;
    const float invDx   = 1.0f / dx;
    const float omega   = 1.9f;                 // over-relax (SOR)

    //SetReflectiveBoundary();
    //SetOpenBoundary();

    int blockH = (rows - 2) / 2;                // rows   per block stripe
    int blockW = (cols - 2) / 4;                // columns per block stripe

    auto divBlock = [&](int r0, int c0, int r1, int c1)
    {
        for (int i = r0; i < r1; ++i) {
            for (int j = c0; j < c1; ++j)
            {
                if (i < 1 || i >= rows-1 || j < 1 || j >= cols-1) continue;

                float du_dx = (uVel.GetValue(i, j+1) - uVel.GetValue(i, j))   * invDx;
                float dv_dy = (vVel.GetValue(i+1, j) - vVel.GetValue(i, j))   * invDx;
                divergence.SetValue(i, j, du_dx + dv_dy);
            }
        }
    };

    {
        std::vector<std::thread> th;
        for (int bi = 0; bi < 2; ++bi)
            for (int bj = 0; bj < 4; ++bj)
            {
                int r0 = 1 +  bi      * blockH;
                int r1 = 1 + (bi + 1) * blockH;
                int c0 = 1 +  bj      * blockW;
                int c1 = 1 + (bj + 1) * blockW;
                th.emplace_back(divBlock, r0, c0, r1, c1);
            }
        for (auto& t : th) t.join();
    }

    pressure = ValueField(rows, cols, 0.0f);
    auto pressureBlock = [&](int r0, int c0, int r1, int c1)
    {
        for (int it = 0; it < gaussSeidelIterations; ++it)
            for (int i = r0; i < r1; ++i)
                for (int j = c0; j < c1; ++j)
                {
                    if (i < 1 || i >= rows-1 || j < 1 || j >= cols-1) continue;

                    float pL  = pressure.GetValue(i,   j-1);
                    float pR  = pressure.GetValue(i,   j+1);
                    float pB  = pressure.GetValue(i-1, j);
                    float pT  = pressure.GetValue(i+1, j);
                    float rhs = divergence.GetValue(i, j);

                    float pNew = 0.25f * (pL + pR + pB + pT - rhs * dx*dx);
                    float pOld = pressure.GetValue(i, j);
                    pressure.SetValue(i, j, pOld + omega * (pNew - pOld));
                }
    };

    {
        std::vector<std::thread> th;
        for (int bi = 0; bi < 2; ++bi)
            for (int bj = 0; bj < 4; ++bj)
            {
                int r0 = 1 +  bi      * blockH;
                int r1 = 1 + (bi + 1) * blockH;
                int c0 = 1 +  bj      * blockW;
                int c1 = 1 + (bj + 1) * blockW;
                th.emplace_back(pressureBlock, r0, c0, r1, c1);
            }
        for (auto& t : th) t.join();
    }

    auto uBlock = [&](int r0, int c0, int r1, int c1)
    {
        for (int i = r0; i < r1; ++i)
            for (int j = c0; j < c1; ++j)
            {
                if (j < 1 || j >= cols) continue;      // u faces skip walls

                float gradP = (pressure.GetValue(i, j) - pressure.GetValue(i, j-1)) * invDx;
                uVel.SetValue(i, j, uVel.GetValue(i, j) - gradP);
            }
    };

    {
        std::vector<std::thread> th;
        for (int bi = 0; bi < 2; ++bi)
            for (int bj = 0; bj < 4; ++bj)
            {
                int r0 =      bi      * blockH;
                int r1 = (bi == 1) ? rows : r0 + blockH;
                int c0 =  bj      * blockW;
                int c1 = (bj == 3) ? cols+1 : c0 + blockW;
                th.emplace_back(uBlock, r0, c0, r1, c1);
            }
        for (auto& t : th) t.join();
    }

    auto vBlock = [&](int r0, int c0, int r1, int c1)
    {
        for (int i = r0; i < r1; ++i)
        {
            if (i < 1 || i >= rows) continue;          // v faces skip walls
            for (int j = c0; j < c1; ++j)
            {
                float gradP = (pressure.GetValue(i, j) - pressure.GetValue(i-1, j)) * invDx;
                vVel.SetValue(i, j, vVel.GetValue(i, j) - gradP);
            }
        }
    };

    {
        std::vector<std::thread> th;
        for (int bi = 0; bi < 2; ++bi)
            for (int bj = 0; bj < 4; ++bj)
            {
                int r0 =  bi      * blockH;
                int r1 = (bi == 1) ? rows+1 : r0 + blockH;
                int c0 = 1 + bj      * blockW;          // same j extents as pressure
                int c1 = 1 + (bj + 1) * blockW;
                th.emplace_back(vBlock, r0, c0, r1, c1);
            }
        for (auto& t : th) t.join();
    }
    //SetReflectiveBoundary();
}

void FluidSolverC::InjectRandomEddiesAtTop(float chance = 0.01f)
{
    for (int j = 1; j < cols - 1; ++j)
    {
        if (rand() / float(RAND_MAX) < chance)
        {
            // Inject downward v-velocity and some curl
            vVel.SetValue(0, j, -3.0f);  // enters downward
            uVel.SetValue(0, j,   1.0f);  // swirl left
            uVel.SetValue(0, j+1, -1.0f); // swirl right

            // Optional: inject dye
            currentDensity.SetValue(1, j, 0.5f);
        }
    }
}

void FluidSolverC::InjectInstantBlobAtTop(float strength = 120.0f, int width = 20, float density = 1.0f)
{
    // Don't inject too close to side walls
    if (cols <= width + 2) return;

    int startJ = rand() % (cols - width - 2) + 1;

    for (int j = startJ; j < startJ + width; ++j)
    {
        // Normalized horizontal position within the blob [0, 1]
        float x = float(j - startJ) / float(width);

        // Bell-curve profile: peak in the center
        float profile = 1.0f - 4.0f * (x - 0.5f) * (x - 0.5f);  // parabola, max = 1

        // === Inject downward velocity (v-component) ===
        vVel.SetValue(1, j, -strength * profile);  // inflow into domain

        // === Swirl: optional horizontal rotation using u-component ===
        float swirl = strength * (x - 0.5f);  // antisymmetric around center
        uVel.SetValue(1, j,   +swirl);
        uVel.SetValue(1, j+1, -swirl);

        // === Optional: inject some density to visualize ===
        currentDensity.SetValue(1, j, density * profile);
    }
}

void FluidSolverC::MaybeSpawnBlob()
{
    static std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<float> rand01(0.0f, 1.0f);
    std::uniform_int_distribution<int> randCol(19, cols - 20);

    if (rand01(rng) < 0.8f) { // 10%
        BlobInflow blob;
        blob.radius   = 3 + int(rand01(rng) * 15);     // 3–6 radius
        blob.strength = 0.6f + rand01(rng) * 0.2f;
        blob.jCenter  = randCol(rng);
        blob.direction = glm::vec2((rand01(rng) - 0.5f) * 0.2f, 1.0f);
        blob.lifeLeft = blob.radius * 2 + 1;
        blob.currentLayer = 0;

        // Precompute circular layers: each row is a 1D slice of the blob
        int diam = blob.radius * 2 + 1;
        blob.layers.resize(diam);

        for (int i = 0; i < diam; ++i) {
            blob.layers[i].resize(diam);
            for (int j = 0; j < diam; ++j) {
                int di = i - blob.radius;
                int dj = j - blob.radius;
                float dist2 = di*di + dj*dj;
                if (dist2 > blob.radius * blob.radius) continue;

                float falloff = 1.0f - dist2 / (blob.radius * blob.radius);
                blob.layers[i][j] = falloff * blob.strength;
            }
        }

        activeBlobs.push_back(blob);
    }
}

// void FluidSolverC::InjectBlob(BlobInflow& blob)
// {
//     int layerIdx = blob.currentLayer;
//     if (layerIdx >= blob.layers.size()) return;

//     const auto& rowSlice = blob.layers[layerIdx];
//     int diam = rowSlice.size();

//     int i = 1;  // inject at top interior row
//     for (int dj = 0; dj < diam; ++dj)
//     {
//         int j = blob.jCenter + dj - blob.radius;
//         if (j < 1 || j >= cols - 1) continue;

//         float value = rowSlice[dj];
//         currentDensity.SetValue(i, j, std::max(currentDensity.GetValue(i, j), value));

//         float velocityBoost = 30.0f;  // or whatever scale makes sense for your system

//         uVel.SetValue(i, j, blob.direction.x * value * velocityBoost);
//         vVel.SetValue(i, j, blob.direction.y * value * velocityBoost);
//     }

//     blob.currentLayer++;
//     blob.lifeLeft--;
// }

void FluidSolverC::InjectBlob(BlobInflow& blob)
{
    int layerIdx = blob.currentLayer;
    if (layerIdx >= blob.layers.size()) return;

    const auto& rowSlice = blob.layers[layerIdx];
    int diam = rowSlice.size();

    int i = 1;  // inject at top interior row
    float velocityBoost = 30.0f;
    float swirlStrength = 0.3f;

    for (int dj = 0; dj < diam; ++dj)
    {
        int j = blob.jCenter + dj - blob.radius;
        if (j < 1 || j >= cols - 1) continue;

        float value = rowSlice[dj];
        if (value < 0.01f) continue;

        currentDensity.SetValue(i, j, std::max(currentDensity.GetValue(i, j), value));

        glm::vec2 mainFlow = blob.direction * value * velocityBoost;

        float offset = dj - blob.radius;
        glm::vec2 swirlDir = glm::normalize(glm::vec2(offset, -1.0f));  // CLOCKWISE

        glm::vec2 swirl = swirlDir * value * swirlStrength * velocityBoost;
        glm::vec2 totalVel = mainFlow + swirl;

        uVel.SetValue(i, j, totalVel.x);
        vVel.SetValue(i, j, totalVel.y);
    }

    blob.currentLayer++;
    blob.lifeLeft--;
}

void FluidSolverC::Step()
{
    float baseVelocity =  60.0f / 1.1f;      // reuse your strength heuristic
    glm::vec2 sourceVel(0.0f, -baseVelocity);

    const int offset = 15;
    int i0 = rows - 1 - offset;
    int i1 = rows - 1;
    int j0 = cols / 2 - offset;
    int j1 = cols / 2 + offset;

    /* Add dye at cell centres */
    for (int i = i0; i < i1; ++i)
        for (int j = j0; j < j1; ++j)
            currentDensity.SetValue(i, j, 1.0f);
    
    int centreCol = cols / 2;
    for (int i = i1; i >= i0; --i)
    {
        uVel.SetValue(i, centreCol, sourceVel.x);   // = 0 in this demo
        uVel.SetValue(i-1, centreCol, sourceVel.x);   // = 0 in this demo

        vVel.SetValue(i, centreCol, sourceVel.y);
        vVel.SetValue(i, centreCol-1, sourceVel.y);
    }

    Diffusion();
    Advection();
    ClearDivergence();

    // if (rand() / float(RAND_MAX) < 0.01f) {
    //     InjectInstantBlobAtTop();
    // }

    //Blob related shenanigans
    const int maxActiveBlobs = 50;
    if (activeBlobs.size() < maxActiveBlobs)
    {
        MaybeSpawnBlob();
    }    

    for (auto it = activeBlobs.begin(); it != activeBlobs.end(); ) {
        InjectBlob(*it); // injects current layer and decreases lifeLeft

        if (it->lifeLeft <= 0 || it->currentLayer >= it->layers.size()) {
            // Blob has finished its life or used up all layers → remove it
            it = activeBlobs.erase(it);
        } else { ++it; }
    }
}

std::vector<float> FluidSolverC::GenerateArrowPositions(float startX,
                                                        float startY,
                                                        float h,
                                                        float maxStrength) const
{
    // interior cell count (one row/col less than the staggered fields)
    const int N = rows - 1;
    const int M = cols - 1;

    std::vector<float> pos;
    pos.reserve(N * M * 6);

    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < M; ++j)
        {
            glm::vec2 v = GetVelocityAtCenter(i, j);
            float     s = glm::length(v) / maxStrength;

            float cx = startX + j * h + 0.5f * h;  // cell centre
            float cy = startY - i * h - 0.5f * h;

            glm::vec2 tip = glm::normalize(v) * 0.5f * h * s;

            pos.insert(pos.end(),
            {
                cx,           cy,           0.0f,   // tail
                cx + tip.x,   cy + tip.y,   0.0f    // head
            });
        }
    }
    return pos;
}

std::vector<float> FluidSolverC::GenerateArrowColours(float maxStrength,
                                                      const glm::vec3& cSlow,
                                                      const glm::vec3& cFast) const
{
    const int N = rows - 1;
    const int M = cols - 1;

    std::vector<float> col;
    col.reserve(N * M * 6);

    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < M; ++j)
        {
            glm::vec2 v = GetVelocityAtCenter(i, j);

            float t = glm::length(v) / maxStrength;
            if (t < 0.0f) t = 0.0f;        // manual clamp
            else if (t > 1.0f) t = 1.0f;

            glm::vec3 c = (1.0f - t) * cSlow + t * cFast;

            // duplicate for tail & head
            for (int k = 0; k < 2; ++k)
                col.insert(col.end(), { c.r, c.g, c.b });
        }
    }
    return col;
}