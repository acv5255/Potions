#pragma once

enum PotionsRunType {
    EQUILIBRIUM = 0,
    KINETIC = 1
};

const double TOLERANCE = 1e-12;  // Tolerance for the numerical solver
const int MAX_ITERS = 50;        // Limit number of iterations