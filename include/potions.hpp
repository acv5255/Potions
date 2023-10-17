#pragma once

#include <stdexcept>
#include "constants.hpp"
#include "potions_input.hpp"
#include "solution.hpp"

class NotImplemented : public std::logic_error {
    public:
        NotImplemented() : std::logic_error("Function not yet implemented") { };
};