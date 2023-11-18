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
    vec totalConcentration;

    ChemicalState() { };
    ChemicalState(const vec& conc, const vec& totConc) : concentration(conc), totalConcentration(totConc) { };
};

struct EquilibriumConstants {
    mat stoichMat;          // Stoichiometry matrix for all equilibrium reactions
    vec eqConsts;           // Equilibrium constants

    EquilibriumConstants() { };
    EquilibriumConstants(const mat& stoich_mat, const vec& eq_consts) : stoichMat(stoich_mat), eqConsts(eq_consts) { };
};

struct TotalConstants {
    mat tot_mat;     // Total concentration matrix

    TotalConstants() { };
    TotalConstants(const mat& tot_mat) : tot_mat(tot_mat) { };
};

struct KineticConstants {
    mat kin_mat;         // Kinetic matrix
    vec kin_const;     // Kinetic rate constants
    vec eq_const;       // Kinetic equilibrium constants

    KineticConstants() { };
    KineticConstants(const mat& kin_mat, const vec& kin_consts, const vec& eq_const) :
        kin_mat(kin_mat),
        kin_const(kin_consts), eq_const(eq_const) { };
};
