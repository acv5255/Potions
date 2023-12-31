#include "ode.hpp"

/*
    Solve a system of first-order ordinary differential equations using the classical 4th-order Runge-Kutta
    method. This method is locally 5th order consistent and 4th order globally consistent. This method
    involves 4 function evalutions and is relatively fast given it's consistency.
 */
template<typename T>
Col<T> runge_kutta_explicit(const function<Col<T>(Col<T>)>& f, const Col<T>& x0, T dt) {
    const double one_sixth = 0.1666666666666666666667;
    const Col<T> k1 = f(x0);
    const Col<T> k2 = f(x0  + 0.5 * dt * k1);
    const Col<T> k3 = f(x0 + 0.5 * dt * k2);
    const Col<T> k4 = f(x0 + dt * k3);

    return x0 + dt * one_sixth * (k1 + k4 + 2.0 * (k2 + k3) );
}

template Col<f64> runge_kutta_explicit<f64>(const function<Col<f64>(Col<f64>)>&, const Col<f64>&, f64);

