#pragma once

#include <stdexcept>
#include "utils.hpp"
#include "constants.hpp"
#include "potions_input.hpp"
#include "solution.hpp"
#include "model_inputs.hpp"

class NotImplemented : public std::logic_error {
    public:
        NotImplemented() : std::logic_error("Function not yet implemented") { };
};