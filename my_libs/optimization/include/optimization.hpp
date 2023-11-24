#pragma once

/** @file */

#include <stdexcept>
#include <functional>
#include <optional>
#include "armadillo"

using std::function;
using std::optional;
using arma::Col;
using arma::Mat;

const unsigned int MAX_ITERATIONS = 50; /** Maximum number of iterations before failing */
using f64 = double;

/**
 * @brief Numerically approximate the jacobian matrix for a vector function
 * Compute the numerical approximation for the Jacobian matrix using the
 * centered-difference approximation. This requires the size of the function
 * 'f' to be the number of input parameters so that the corresponding Jacobian
 * matrix is square.
 * @tparam T A floating point data type
 * @param func The continuous function whose Jacobian to approximate
 * @param x The value at which to approximate the Jacobian of 'f'
 * @return Mat<T> 
 */
template<typename T>
Mat<T> jacobian(const function<Col<T>(Col<T>)>& func, const Col<T>& x);

/**
 * @brief Determine the root of a general vector function
 * Determine the root of a general vector function using Newton-Raphson iteration.
 * Newton-Raphson iteration requires evaluating a Jacobian and then solving the 
 * resulting linear system. The iterative method looks like:
 * x^{i+1} = x^{i} - J(x^{i})^{-1}f(x^{i})
 * This method has quadratic convergence, though may be unstable and therefore requires
 * an initial guess that is sufficiently close to the root.
 * @tparam T 
 * @param func The function to find the root of
 * @param jac The function that generates the jacobian matrix
 * @param x0 The initial guess for the root
 * @return optional<Col<T>> 
 */
template<typename T>
optional<Col<T>> root(const function<Col<T>(Col<T>)>& func, const function<Mat<T>(Col<T>)>& jac, const Col<T>& x0);