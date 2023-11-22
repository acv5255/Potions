#pragma once

#include <string>
#include "constants.hpp"

using std::string;

class RunOptions {
    string _output_name;          // Name of the output file to write

    public:
        RunOptions() { };
        string output_name() {return _output_name; }

    static RunOptions from_arguments(int argc, char* argv[]);
};

