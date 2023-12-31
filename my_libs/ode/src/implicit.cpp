#include "ode.hpp"

/*
    Solves the differential equation system described by `f` using the 4th-order consistent
    Runge-Kutta-Gauss method, an implicit collocation method. Because it is an implicit method,
    this solver requires numerical optimization, which is implemented by my other package `optimization`
    using Newton-Raphson iteration.
 */
template<typename T>
Col<T> runge_kutta_implicit(const function<Col<T>(Col<T>)>& f, const Col<T>& x0, T dt) {
    /*
        Solve the autonomous system using Runge-Kutta method
     */
    const int N = x0.size();

    // Define the Runge-Kutta-Gauss parameters
    const double b1 = 0.5;
    const double b2 = 0.5;

    const double a11 = 0.25;
    const double a12 = 0.25 - std::sqrt(3) / 6.0;
    const double a21 = 0.25 + std::sqrt(3) / 6.0;
    const double a22 = 0.25;

    // Construct the residual functions for rootfinding
    auto k1_res = [&] (const Col<T>& k1, const Col<T>& k2)  {
        Col<T> val = f(x0 + a11 * k1 * dt + a12 * k2 * dt) - k1;
        return val;
    };

    auto k2_res = [&] (const Col<T>& k1, const Col<T>& k2) {
        Col<T> val = f(x0 + a21 * k1 * dt + a22 * k2 * dt) - k2;
        return val;
    };

    Col<T> ks_0 = arma::zeros(2 * N);
    Col<T> k_guess = f(x0);
    for (int i = 0; i < N; i++) {
        ks_0[i] = k_guess[i];
        ks_0[i + N] = k_guess[i];
    }

    // Combine the two residual vectors into a single function to provide the numerical optimizer
    function<Col<T>(Col<T>)> residual = [&] (const Col<T>& ks) -> Col<T> {
        Col<T> res = arma::zeros(2 * N);
        Col<T> k1 = arma::zeros(N);
        Col<T> k2 = arma::zeros(N);
        for (int i = 0; i < N; i++) {
            k1[i] = ks[i];
            k2[i] = ks[i+N];
        }

        Col<T> k1_val = k1_res(k1, k2);
        Col<T> k2_val = k2_res(k1, k2);

        for (int i = 0; i < N; i++) {
            res[i] = k1_val[i];
            res[i + N] = k2_val[i];
        }

        return res;
    };

    // Define the jacobian function for the Newton-Raphson iteration
    auto jac = [&] (const Col<T>& val) {
        return jacobian<T>(residual, val);
    };

    // Need an initial gues for the solution
    Col<T> res_k0 = residual(ks_0);

    // Now, find the root of this function
    optional<Col<T>> ks_root_option = root<T>(residual, jac, ks_0);

    if (!ks_root_option.has_value()) {
        std::cerr << "Error: Failed to find root in Implicit Runge-Kutta solver\n";
        exit(-1);
    }

    Col<T> ks = ks_root_option.value();

    Col<T> k1 = arma::zeros(N);
    Col<T> k2 = arma::zeros(N);

    for (int i = 0; i < N; i++) {
        k1[i] = ks[i];
        k2[i] = ks[i + N];
    }

    // Determine the next step with the calculated Runge-Kutta parameters
    return x0 + dt * (b1 * k1 + b2 * k2);
}

template Col<f64> runge_kutta_implicit<f64>(const function<Col<f64>(Col<f64>)>&, const Col<f64>&, f64);