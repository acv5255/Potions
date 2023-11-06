#pragma once
/*
    Functions for solving equilibrium
 */
#include "armadillo"
#include "potions.hpp"

template<typename T>
ChemicalState<T> SolveEquilibrium(const Col<T>& tot_conc, const EquilibriumConstants<T>& eq_params, const TotalConstants<T>& tot_params);