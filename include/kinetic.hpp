#pragma once

#include "potions.hpp"

/*
    Return the kinetic rates for each species in moles per second
 */
vec kinetic_rate(const vec& conc, const vec& surface_area, const KineticConstants& kin_params);

// template<typename T>
