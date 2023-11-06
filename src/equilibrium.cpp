#include <functional>
#include <iostream>
#include "potions.hpp"
#include "optimization.hpp"

using std::function;

template<typename T>
ChemicalState<T> SolveEquilibrium(const Col<T>& tot_conc, const EquilibriumConstants<T>& eq_params, const TotalConstants<T>& tot_params) {
    const Col<T> log_x_p = arma::pinv(eq_params.stoich_mat) * eq_params.eq_consts;  // Particular solution
    const Mat<T> stoich_mat_nullspace = arma::null(eq_params.stoich_mat);
    
    auto log_c = [&] (const Col<T>& log_x) {
        return log_x_p + stoich_mat_nullspace * log_x;
    };

    auto mass_err = [&] (const Col<T>& log_x) {
        const Col<T> c = arma::exp(log_c(log_x));
        return tot_params.tot_mat * c;
    };

    auto jac = [&] (const Col<T> x) -> Mat<T> {
        /*
            Jacobian function
         */
        return jacobian(mass_err, x);
    };

    // Do the optimization now
    const optional<Col<T>> sln = root(mass_err, jac, tot_conc);

    if (sln.has_value) {
        return sln.value;
    }
    else {
        std::cerr << "Equilibrium solution failed to converge\n";
        exit(-1);
    }
}