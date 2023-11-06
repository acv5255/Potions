#include "ode.hpp"

template<typename T>
Col<T> RungeKuttaExplicit(const function<Col<T>(Col<T>)>& func, const Col<T> x0, double dt) {
    const double one_sixth = 0.1666666666666666666667;
    const Col<T> k1 = f(x0);
    const Col<T> k2 = f(x0  + 0.5 * dt * k1);
    const Col<T> k3 = f(x0 + 0.5 * dt * k2);
    const Col<T> k4 = f(x0 + dt * k3);

    return x0 + dt * one_sixth * (k1 + k4 + 2.0 * (k2 + k3) );
}

