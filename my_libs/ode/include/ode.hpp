#include <functional>
#include <stdexcept>
#include "armadillo"
#include "optimization.hpp"

class NotImplemented : public std::logic_error {
    public:
        NotImplemented() : std::logic_error("Function not yet implemented") { };
};

using arma::Col;
using std::function;

const double TOLERANCE = 1e-10;  // Tolerance for implicit ODE


// Solve the ODE using the 4th order method
template<typename T>
Col<T> RungeKuttaExplicit(const function<Col<T>(Col<T>)>& func, const Col<T> x0, double dt);


// Solve an ODE using an implicit Runge-Kutta method
template<typename T>
Col<T> RungeKuttaImplicit(const function<Col<T>(Col<T>)>& func, const Col<T> x0, double dt);

// Solve an ODE using either explicit or implicit method
template<typename T>
Col<T> SolveODE(const function<Col<T>(Col<T>)>& func, const Col<T> x0, double dt);
