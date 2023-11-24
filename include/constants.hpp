#pragma once

/** @file */

/**
 * @brief Type of simulation to run
 * 
 */
enum PotionsRunType {
    EQUILIBRIUM = 0, /* Determine equilibrium at a single time point */
    KINETIC = 1 /* Solve concentrations and multiple time steps with kinetic reactions */
};

const double TOLERANCE = 1e-10; /* Absolute error tolerance for the numerical solvers */
const int MAX_ITERS = 50;       /* Maximum number of iterations for numerical solvers */