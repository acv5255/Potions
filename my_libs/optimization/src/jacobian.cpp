#include "optimization.hpp"

const double DX = 0.01;

template<typename T>
Mat<T> jacobian(const function<Col<T>(Col<T>)>& func, const Col<T>& x) {
    /*
        Approximate the Jacobian matrix for the function `func` at 
        the vector `x`
     */
    const int n = x.size();

    Mat<T> jac = arma::zeros(n, n);

    auto copyVec = [&] () -> Col<T> {
        Col<T> c = Col<T>(x.size());

        for (int i = 0; i < x.size(); i++) {
            c[i] = x[i];
        }

        return c;
    };

    // Compute the jacobian element-wise
    for (int i = 0; i < n; i++) {
        // double dx_i = std::abs(0.01 * x[i]);

        Col<T> x_up = copyVec();
        Col<T> x_dn = copyVec();

        x_up[i] += DX;
        x_dn[i] -= DX;

        Col<T> f_up = func(x_up);
        Col<T> f_dn = func(x_dn);

        Col<T> deriv_x = (f_up - f_dn) / (2 * DX);

        jac.col(i) = deriv_x;

    }

    return jac;
}

template Mat<f64> jacobian<f64>(const function<Col<f64>(Col<f64>)>&, const Col<f64>&);