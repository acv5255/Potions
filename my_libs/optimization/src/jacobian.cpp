#include "optimization.hpp"

template<typename T>
Mat<T> jacobian(const function<Col<T>(Col<T>)>& func, const Col<T>& x) {
    /*
        Approximate the Jacobian matrix for the function `func` at 
        the vector `x`
     */
    const int n = x.size();

    Mat<T> jac = arma::zeros(n, n);

    // Compute the jacobian element-wise

    for (int i = 0; i < n; i++) {
        double dx_i = 0.01 * x[i];
        Col<T> x_up = Col(x);
        Col<T> x_dn = Col(x);
        x_up[i] += dx_i;
        x_dn[i] -= dx_i;

        Col<T> deriv_x = (0.5 / dx_i) * (x_up - x_dn);

        jac.col(i) = deriv_x;

    }

    return jac;
}