#include <cmath>
#include "ode.hpp"

using arma::cx_vec;
using arma::cx_double;

template<typename T>
Col<T> solve_ode(const function<Col<T>(Col<T>)>& func, const Col<T>& x0, T dt) {
    /*
        Solve the ODE using a one-step method. This method employs either either the 
        fourth order explicit Runge-Kutta 4th order method or the 4th order implicit
        Runge-Kutta-Gauss method to solve the system.

        We require the eigenvalues of the Jacobian matrix of the kinetic rates
        to determine the stability of the system.
        We want to use the explicit Runge-Kutta method because
        it is much faster to evaluate, however it cannot solve stiff systems.
        If the time step is not stable, then use the slower, but unconditionally
        stable implicit rule. 
     */
    const Mat<T> jac = jacobian<T>(func, x0);   // <- this is the problem line
    const cx_vec eigen_values = arma::eig_gen(jac);
    
    // Need to check each eigenvalue for the stability
    bool use_explicit = true;
    for (cx_double e_i: eigen_values) {
        const cx_double z = e_i * dt;
        const cx_double rat = 1.0 + z + 0.5 * z * z + (1.0/6.0) * z * z * z + (1.0/24.0) * z * z  * z * z;
        const double stability = std::sqrt(
            rat.real() * rat.real() + rat.imag() * rat.imag()
        );

        /*
            The stability condition of any Runge-Kutta methods 
            requires that this value is less than 1 or else the 
            solution will be unstable. In that case, we must use 
            an unconditionally stable method like the Implicit 
            Runge-Kutta-Gauss method
         */
        if (stability > 1.0) {
            use_explicit = false;
            break;
        }
    }

    if (use_explicit) {
        // std::cerr << "Using explicit method\n";
        return runge_kutta_explicit(func, x0, dt);
    }
    else {
        // std::cerr << "Using implicit method\n";
        return runge_kutta_implicit(func, x0, dt);
    }
}

template Col<f64> solve_ode<f64>(const function<Col<f64>(Col<f64>)>& func, const Col<f64>& x0, f64 dt);