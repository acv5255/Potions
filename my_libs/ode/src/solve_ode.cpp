#include <cmath>
#include "ode.hpp"

using arma::cx_vec;

template<typename T>
Col<T> SolveODE(const function<Col<T>(Col<T>)>& func, const Col<T> x0, double dt) {
    /*
        Solve the ODE by stiffness detection. If the time step is not stable, then
        use the implicit rule
     */
    const Mat<T> jac = jacobian<T>(func, x0);
    const cx_vec eigen_values = arma::eig_gen(jac);
    
    // Get the maximum modulus of the eigenvalue
    bool use_explicit = true;
    for (auto e_i: eigen_values) {
        const double mod = std::sqrt(
            e_i.imag * e_i.imag + e_i.real * e_i.real
        );

        const double z = mod * dt;

        // Stability condition for the 4th order Runge-Kutta method
        const double stability = 
            1.0 + z + 0.5 * z * z + (1.0/6.0) * std::pow(z,3.0) + (1.0/24.0) * std::pow(z,4.0);

        if (stability > 1.0) {
            use_explicit = false;
            break;
        }
    }

    if (use_explicit) {
        return RungeKuttaExplicit(func, x0, dt);
    }
    else {
        return RungeKuttaImplicit(func, x0, dt);
    }
}