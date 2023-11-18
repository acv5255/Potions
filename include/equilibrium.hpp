#pragma once
/*
    Functions for solving equilibrium
 */
#include "armadillo"
#include "potions.hpp"

ChemicalState SolveEquilibrium(const vec& tot_conc, const EquilibriumConstants& eq_params, const TotalConstants& tot_params);