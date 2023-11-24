#pragma once
/** @file */
#include "potions.hpp"

/**
 * @brief Rate of change of species in terms of moles per liter
 * Use the Transition-State-Theory type rate law to determine the rate
 * of change of concentrations of each species in the system. This includes
 * both primary and secondary species. Note: each evaluation of this functions
 * involves solving equilibrium.
 * @param tot_conc The vector of total concentrations with the length being the number of species
 * @param surface_area The vector of mineral surface areas with the length being the number of minerals
 * @param eq The equilibrium parameters for this system
 * @param kin_params The kinetic matrix, rate constants, and kinetic equilibrium constants
 * @param tot_params The mass conservation matrix wrapper structure
 * @return vec The rate of changes of species with the length being the number of species.
 */
vec kinetic_rate(const vec& tot_conc, const vec& surface_area, const EquilibriumConstants& eq, const KineticConstants& kin_params, const TotalConstants& tot_params);

/** @fn
 * @brief Determine time derivative of total concentrations
 * Determine the equilibrium concentrations at time step 'dt' after
 * solving the ODE. This function is computationally expensive and requires
 * solving the autonomous ODE, which will involve multiple evaluations of 
 * 'solve_equilibrium'.
 * 
 * @param chem The chemical state object at the start of the time step
 * @param surfaceArea The vector of mineral surface areas at the time step
 * @param kin The kinetic parameters for our system
 * @param eq The equilibrium parameters for our system
 * @param tot The mass conservations parameters
 * @param dt The final time step to evaluate the ODE
 * @return ChemicalState The chemical state at time dt in the future
 */
ChemicalState solve_kinetic_equilibrium(
    const ChemicalState& chem,
    const vec& surfaceArea,
    const KineticConstants& kin,
    const EquilibriumConstants& eq,
    const TotalConstants& tot,
    double dt
);