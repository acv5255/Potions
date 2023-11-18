#include "potions.hpp"

RunOptions RunOptions::FromArguments(int argc, char* argv[]) {
    RunOptions options;
    if (argc == 2) {
        options._outputName = string(argv[1]);
    }
    else {
        std::cerr << "Expected only input name, other functions not implemented\n";
        exit(-1);
    }

    options._writeOutputs = true;
    options._plotKineticOutputs = true;

    return options;
}

string RunOptions::outputName() {
    return _outputName;
}

bool RunOptions::writeOutputs() {
    return _writeOutputs;
}

bool RunOptions::plotKineticOutputs() {
    return _plotKineticOutputs;
}