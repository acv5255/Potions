#include "potions.hpp"

RunOptions RunOptions::from_arguments(int argc, char* argv[]) {
    RunOptions options;
    if (argc == 2) {
        options._output_name = string(argv[1]);
    }
    else {
        std::cerr << "Expected only input name, other functions not implemented\n";
        exit(-1);
    }

    options._write_outputs = true;
    options._plot_kinetic_outputs = true;

    return options;
}

string RunOptions::output_name() {
    return _output_name;
}

bool RunOptions::write_outputs() {
    return _write_outputs;
}

bool RunOptions::plot_kinetic_outputs() {
    return _plot_kinetic_outputs;
}