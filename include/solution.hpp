#pragma once

#include <tuple>
#include "armadillo"
#include "model_inputs.hpp"

using std::tuple;
using arma::Col;
using arma::Mat;

template<typename T>
struct ChemicalState {
    Col<T> concentration;
    Col<T> totalConcentration;

    ChemicalState(const Col<T>& conc, Col<T>& totConc) : concentration(conc), totalConcentration(totConc) { };
};

template<typename T>
struct EquilibriumConstants {
    Mat<T> stoich_mat;          // Stoichiometry matrix for all equilibrium reactions
    Col<T> eq_consts;           // Equilibrium constants

    EquilibriumConstants(const Mat<T>& stoich_mat, const Col<T>& eq_consts) : stoich_mat(stoich_mat), eq_consts(eq_consts) { };
};

template<typename T>
struct TotalConstants {
    Mat<T> tot_mat;     // Total concentration matrix

    TotalConstants(const Mat<T>& tot_mat) : tot_mat(tot_mat) { };
};

template<typename T>
struct KineticConstants {
    Mat<T> kin_mat;         // Kinetic matrix
    Col<T> rate_consts;     // Kinetic rate constants
    Col<T> eq_consts;       // Kinetic equilibrium constants

    KineticConstants(const Mat<T>& kin_mat, const Col<T>& rate_consts) :
        kin_mat(kin_mat),
        rate_consts(rate_consts) { };
};
