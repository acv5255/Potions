#include <functional>
#include <iostream>
#include "potions.hpp"
#include "optimization.hpp"

using std::function;

ChemicalState solve_equilibrium(const vec& tot_conc, const EquilibriumConstants& eq_params, const TotalConstants& tot_params) {
    // const vec log_x_p = arma::pinv(eq_params.stoichiometry_matrix) * eq_params.equilibrium_constants;  // Particular solution

    auto log_c = [&] (const vec& log_x) -> vec {
        return eq_params.conc_particular_solution + eq_params.stoich_null_space * log_x;
    };

    auto mass_err = [&] (const vec& log_x) -> const vec {
        // const vec c = arma::pow(log_c(log_x), 10.0);
        const vec c = arma::exp(log_c(log_x));
        return tot_params.tot_mat * c - tot_conc;
    };

    auto jac = [&] (const vec& x) -> mat {
        /*
            Jacobian function for the concentration residual
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