#pragma once

#include "potions.hpp"

/*
    Return the kinetic rates for each species in moles per second
 */
vec kinetic_rate(const vec& totConc, const vec& surface_area, const EquilibriumConstants& eq, const KineticConstants& kin_params, const TotalConstants& tot_params);

ChemicalState solve_kinetic_equilibrium(
    const ChemicalState& chem,
    const vec& surfaceArea,
    const KineticConstants& kin,
    const EquilibriumConstants& eq,
    const TotalConstants& tot,
    double dt
);