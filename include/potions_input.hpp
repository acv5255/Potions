#pragma once

#include <string>
#include "constants.hpp"

using std::string;

class RunOptions {
    string _output_name;          // Name of the output file to write
    bool _write_outputs;          // Whether to write the outputs to a file or not
    bool _plot_kinetic_outputs;    // Whether or not to plot the kinetic outputs

    public:
        RunOptions() { };
        string output_name();
        bool write_outputs();
        bool plot_kinetic_outputs();

    static RunOptions from_arguments(int argc, char* argv[]);
};

