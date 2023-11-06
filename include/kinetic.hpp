#pragma once

#include "potions.hpp"

/*
    Return the kinetic rates for each species in moles per second
 */
template<typename T>
Col<T> kinetic_rate(const Col<T>& conc, const Col<T> surface_area, const KineticConstants<T>& kin_params);

// template<typename T>
