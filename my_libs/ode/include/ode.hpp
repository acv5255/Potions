#pragma once

/** @file */

/*
    Author: Andrew Vierbicher
    This file contains the functions used for numerically solving 
    systems of first-order ordinary differential equations.
 */

#include <functional>
#include "armadillo"
#include "optimization.hpp"

using arma::Col;
using std::function;

/**
 * @brief Evolve an ODE using a 4th order explicit Runge-Kutta method
 * Solve a first-order system of ordinary differential equations using the classical
 * 4th-order Runge-Kutta one-step method. Note: this method may not be suitable for stiff
 * systems, though it is much faster than implicit methods.
 * @tparam T 
 * @param func The differential equations to evaulate
 * @param x0 The initial state of the system
 * @param dt The time step to simulate the model
 * @return Col<T> 
 */
template<typename T>
Col<T> runge_kutta_explicit(const function<Col<T>(Col<T>)>& func, const Col<T>& x0, T dt);


/**
 * @brief Evolve an ODE system using the implicit 4th order Runge-Kutta-Gauss method
 * Solve a system of ordinary differential equations using the 4th order Runge-Kutta-Gauss 
 * method. This method is unconditionally stable and will yield a stable solution for stable
 * ODE's, though it is an implicit method and will be slower than explicit methods.
 * @tparam T 
 * @param func The ODE system
 * @param x0 The initial state of the system
 * @param dt The time step
 * @return Col<T> 
 */
template<typename T>
Col<T> runge_kutta_implicit(const function<Col<T>(Col<T>)>& func, const Col<T>& x0, T dt);

/**
 * @brief Solve an ODE system using an explicit method if stable and an unconditionally stable implicit Runge-Kutta method instead.
 * Solve an ODE over a single time step using either an implicit or explicit method without
 * having to worry about the stiffness. Both the implicit and explicit methods are 4th order
 * consistent. This method will determine the stability of the ODE to determine the best solver. 
 * If the system is stable, the method will solve with the faster Runge-Kutta method, and will
 * use the slower, yet stable, implicit method if the explicit method is not stable.
 * @tparam T 
 * @param func The ODE system
 * @param x0 The initial state of the system
 * @param dt The time step
 * @return Col<T> 
 */
template<typename T>
Col<T> solve_ode(const function<Col<T>(Col<T>)>& func, const Col<T>& x0, T dt);
