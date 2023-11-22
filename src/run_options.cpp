#include "potions.hpp"

/*
    This function reads the command line parameters and constructs
    the run options structure. This is a bare-bones data structure,
    but it is set up so that it can be edited in the future.
 */
RunOptions RunOptions::from_arguments(int argc, char* argv[]) {
    RunOptions options;
    if (argc == 2) {
        options._output_name = string(argv[1]);
    }
    else {
        std::cerr << "Expected only input name, other functions not implemented\n";
        exit(-1);
    }

    return options;
}
