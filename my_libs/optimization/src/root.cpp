#include "optimization.hpp"

template<typename T>
optional<Col<T>> root(const function<Col<T>(Col<T>)>& func, const function<Mat<T>(Col<T>)>& jac, const Col<T>& x0) {
    /*
        Determine the root of the function
     */
    Col<T> x = Col<T>(x0);

    function<T(Col<T>)> err = [&] (const Col<T>& x) -> T {
        Col<T> x_sqr = arma::pow(func(x), 2.0);
        T mean_val = arma::sum(x_sqr) / x.size();
        return mean_val;
    };

    for (int i = 0; i < MAX_ITERATIONS; i++) {
        const Col<T> dx = arma::solve(jac(x), func(x));
        x = x + (-dx);

        if (err(x) < TOLERANCE) {
            return optional<Col<T>>(x);
        }
    }

    return std::nullopt;
}