#pragma once

#include <string>
#include "constants.hpp"

using std::string;

class RunOptions {
    string _outputName;          // Name of the output file to write
    bool _writeOutputs;          // Whether to write the outputs to a file or not
    bool _plotKineticOutputs;    // Whether or not to plot the kinetic outputs

    public:
        string outputName();
        bool writeOutputs();
        bool plotKineticOutputs();

    static RunOptions FromArguments(int argc, char* argv[]);
};

