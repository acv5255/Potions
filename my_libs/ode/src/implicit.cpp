#include <stdexcept>
#include "ode.hpp"

class NotImplemented : public std::logic_error {
    public:
        NotImplemented() : std::logic_error("Function not yet implemented") { };
};

template<typename T>
Col<T> RungeKuttaImplicit(const function<Col<T>(Col<T>)>& func, const Col<T> x0, T dt) {
    /*
        Solve the 
     */
    // const double b1 = 

    throw NotImplemented();
}