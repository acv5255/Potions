#pragma once

#include <stdexcept>
#include "constants.hpp"
#include "potions_input.hpp"
#include "solution.hpp"
#include "model_inputs.hpp"
#include "kinetic.hpp"
#include "equilibrium.hpp"
#include "potions_utils.hpp"

class NotImplemented : public std::logic_error {
    public:
        NotImplemented() : std::logic_error("Function not yet implemented") { };
};