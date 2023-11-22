#include <iostream>
#include "optimization.hpp"

const double TOLERANCE = 1e-7;

template<typename T>
void print_matrix(const Mat<T>& m) {
    for (int i = 0; i < m.n_rows; i++) {
        for (int j = 0; j < m.n_cols; j++) {
            std::cout << m(i,j) << " ";
        }

        std::cout << "\n";
    }
}

template<typename T>
optional<Col<T>> root(const function<Col<T>(Col<T>)>& func, const function<Mat<T>(Col<T>)>& jac, const Col<T>& x0) {
    /*
        Determine the root of the function
     */
    Col<T> x = arma::zeros(x0.size());
    for (int i = 0; i < x.size(); i++) x[i] = x0[i];

    function<T(Col<T>)> err = [&] (const Col<T>& val) -> T {
        Col<T> x_sqr = arma::pow(func(val), 2.0);
        T mean_val = arma::sum(x_sqr) / val.size();
        return mean_val;
    };

    double max_err = 9999.0;

    for (int i = 0; i < MAX_ITERATIONS; i++) {
        Mat<T> jac_mat = jac(x);

        const Col<T> dx = arma::solve(jac(x), func(x));
        x = x + (-dx);

        max_err = err(x);

        if (max_err < TOLERANCE) {
            return optional<Col<T>>(x);
        }
    }
    std::cout << "Reached maximum iterations\n";

    return std::nullopt;
}

template optional<Col<f64>> root<f64>(const function<Col<f64>(Col<f64>)>& func, const function<Mat<f64>(Col<f64>)>& jac, const Col<f64>& x0);