#pragma once

#include <tuple>
#include "armadillo"

using std::tuple;
using arma::Col;
using arma::Mat;
using arma::vec;
using arma::mat;

struct ChemicalState {
    vec concentration;
    vec total_concentration;

    ChemicalState() { }
    ChemicalState(const vec& conc, const vec& tot_conc) : concentration(conc), total_concentration(tot_conc) { }
};

struct EquilibriumConstants {
    mat stoichiometry_matrix;          // Stoichiometry matrix for all equilibrium reactions
    mat stoich_null_space;      // Null space used for equilibrium solution
    vec conc_particular_solution;
    vec equilibrium_constants;           // Equilibrium constants

    EquilibriumConstants() { }
    EquilibriumConstants(const mat& stoich_mat, const vec& eq_consts);
};

struct TotalConstants {
    mat tot_mat;     // Total concentration matrix

    TotalConstants() { }
    TotalConstants(const mat& tot_mat) : tot_mat(tot_mat) { };
};

struct KineticConstants {
    mat kin_mat;       // Kinetic matrix
    vec kin_const;     // Kinetic rate constants
    vec eq_const;      // Kinetic equilibrium constants

    KineticConstants() { };
    KineticConstants(const mat& kin_mat, const vec& kin_consts, const vec& eq_const) :
        kin_mat(kin_mat),
        kin_const(kin_consts),
        eq_const(eq_const) { }
};
