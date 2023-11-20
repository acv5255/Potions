#include <iostream>
#include "optimization.hpp"

const double TOLERANCE = 1e-7;

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

    double max_err = 9999.0;

    for (int i = 0; i < MAX_ITERATIONS; i++) {
        // std::cout << "Starting iteration " << i << "\n";
        Mat<T> jac_mat = jac(x);

        // std::cout << "Jacobian: \n";
        // for (int i = 0; i < jac_mat.n_rows; i++) {
        //     for (int j = 0; j < jac_mat.n_cols; j++) {
        //         std::cout << jac_mat(i,j) << " ";
        //     }
        //     std::cout << "\n";
        // }

        const Col<T> dx = arma::solve(jac(x), func(x));
        x = x + (-dx);

        max_err = err(x);
        std::cout << "Error after iteration " << i << ": " << max_err << "\n";

        if (max_err < TOLERANCE) {
            return optional<Col<T>>(x);
        }
    }

    return std::nullopt;
}

template optional<Col<f64>> root<f64>(const function<Col<f64>(Col<f64>)>& func, const function<Mat<f64>(Col<f64>)>& jac, const Col<f64>& x0);