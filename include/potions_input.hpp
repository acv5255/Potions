#pragma once

/** @file */

#include <string>
#include "constants.hpp"

using std::string;

/**
 * @brief Wrapper class for the command line options
 */
class RunOptions {
    string _output_name;    /** Name of the input folder name */

    public:
        /**
         * @brief Construct a new Run Options object
         * 
         */
        RunOptions() { }

        /**
         * @brief Getter for the input folder name
         * 
         * @return string 
         */
        string output_name() {return _output_name; }

        /**
         * @brief Directly handle the command line parameters
         * 
         * @param argc The number of command line arguments
         * @param argv Arguments passed to the command line
         * @return RunOptions 
         */
        static RunOptions from_arguments(int argc, char* argv[]);
};

