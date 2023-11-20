#include <functional>
#include <iostream>
#include "potions.hpp"
#include "optimization.hpp"

using std::function;

ChemicalState SolveEquilibrium(const vec& tot_conc, const EquilibriumConstants& eq_params, const TotalConstants& tot_params) {
    const vec log_x_p = arma::pinv(eq_params.stoichMat) * eq_params.eqConsts;  // Particular solution
    
    // std::cout << "Particular solution: \n";
    // for (auto x: log_x_p) std::cout << x << std::endl;
    // std::cout << std::endl;

    const mat stoich_mat_nullspace = arma::null(eq_params.stoichMat);

    // std::cout << "Null space: \n";
    // for (int i = 0; i < stoich_mat_nullspace.n_rows; i++) {
    //     for (int j = 0; j < stoich_mat_nullspace.n_cols; j++) {
    //         std::cout << stoich_mat_nullspace(i,j) << " ";
    //     }
    //     std::cout << "\n";
    // }
    
    // std::cout << std::endl;

    auto log_c = [&] (const vec& log_x) -> vec {
        return log_x_p + stoich_mat_nullspace * log_x;
    };

    auto mass_err = [&] (const vec& log_x) -> const vec {
        // const vec c = arma::pow(log_c(log_x), 10.0);
        const vec c = arma::exp(log_c(log_x));
        return tot_params.tot_mat * c - tot_conc;
    };

    auto jac = [&] (const vec& x) -> mat {
        /*
            Jacobian function
         */
        return jacobian<double>(mass_err, x);
    };



    // Do the optimization now
    vec xInit = arma::zeros(tot_conc.size());
    const optional<vec> sln = root<double>(mass_err, jac, tot_conc);

    if (sln.has_value()) {
        // vec concVec = arma::pow(log_c(sln.value()), 10.0);
        vec concVec = arma::exp(log_c(sln.value()));

        return ChemicalState(
            concVec,
            tot_conc
        );
    }
    else {
        std::cerr << "Equilibrium solution failed to converge\n";
        exit(-1);
    }
}