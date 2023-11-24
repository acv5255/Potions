#pragma once
#include "armadillo"
#include "potions.hpp"

/** @file */

/** 
 * @brief Solve for the thermodynamic equilibrium
 * This functions solves for the equilibrium speciation for primary and secondary species.
 * In general solving for equilibrium is a highly nonlinear process that involves numerical
 * optimization. The number of equations that we must solve is equal to the number of primary
 * species plus the number of secondary species, and the governing equations come from
 * satisfying both a mass balance and the equilibrium relations:
 * - For ever primary species, there is one mass conservation relationship*
 *  - *For H+, we use a charge balance rather than a conservation of H+
 * - For every secondary species, there is an equilibrium relationship to satisfy
 * 
 * @param tot_conc The vector of total concentrations for each species for the mass balance.
 * @param eq_params The stoichiometry matrix and equilibrium constants for the chemical system
 * @param tot_params The mass concentration matrix
 * @return ChemicalState The ChemicalState object of equilibrium concentrations and a copy of the total concentrations at equilibrium.
 */
ChemicalState solve_equilibrium(const vec& tot_conc, const EquilibriumConstants& eq_params, const TotalConstants& tot_params);