#include "potions.hpp"

EquilibriumConstants::EquilibriumConstants(const mat& stoich_mat, const vec& eq_consts) {
    this->stoichiometry_matrix = stoich_mat;
    this->stoich_null_space = arma::null(stoich_mat);
    this->conc_particular_solution = arma::pinv(stoich_mat) * eq_consts;
    this->equilibrium_constants = eq_consts;
}