#include <vector>
#include <string>
#include <filesystem>
#include <utility>
#include "model_inputs.hpp"

/** @file */

using std::string;
using std::filesystem::path;
using std::vector;
using std::pair;

/**
 * @brief Construct the output file path
 * Construct the output file path to write the simulation results to. The file path
 * has the format "<sim_name>_YYYYMMDD_HHMMSS_results.csv", where YYYYMMDD_HHMMSS is
 * a time stamp for the simulation time and 'sim_name' is the name of the input folder.
 * This is intended so that the results name will always be unique and always be
 * idenfitiable. The results are written in a comma-separated format, so the file is
 * a CSV file.
 * @param simulation_name The name of the input folders
 * @return string 
 */
string get_output_file_path(string simulation_name);

/**
 * @brief Compare floating point numbers safely
 * 
 * @param a 
 * @param b 
 * @return true 
 * @return false 
 */
bool compare_doubles(double a, double b);

/**
 * @brief Write the results of an equilibrium simulation to a file
 * This function writes the equilibrium concentrations to a CSV file with
 * 2 columns:
 *  1. The name of each of the species
 *  2. The concentration (in molars) of each species
 * @param chms 
 * @param species 
 * @param file_path 
 * @return true 
 * @return false 
 */
bool save_equilibrium_results(const ChemicalState& chms, vector<string> species, const string& file_path);

/**
 * @brief Write the results of a kinetic simulation to a CSV file
 * This function writes the results of the simulation to a CSV file as a table.
 * Only the concentrations, not total concentrations, are written to the file.
 * The first row of the CSV file are the column names. The first is the time (in seconds),
 * followed by each chemical species name.
 * The columns in this file are:
 *  1. The time (in seconds) for the simulation step
 *  2. From now on, the concentrations of each species
 * @param res 
 * @param species 
 * @param filePath 
 * @return true 
 * @return false 
 */
bool save_kinetic_results(const vector<pair<double, ChemicalState>>& res, const vector<string>& species, const string& filePath);

/**
 * @brief Get the ionic charge of a chemical species
 * Determine the ionic charge of a chemical species by its formula.
 * Example:
 *  get_charge("Cl-") -> -1
 *  get_charge("Ca++") -> 2
 *  get_charge("H2CO3") -> 0
 * @param name 
 * @return int 
 */
int  get_charge(const string& name);

/**
 * @brief Print the contents of a matrix to the standard output
 * 
 * @tparam T 
 * @param m 
 */
template<typename T>
void print_matrix(const Mat<T>& m);

/**
 * @brief Produce a plot of the kinetic simulation results
 * Produce a time series line plot of the concentrations of each chemical species
 * from a kinetic simulation.
 * @param results The vector of pairs of simulation step results
 * @param species_names The vector of species names
 * @return true 
 * @return false 
 */
bool plot_results(const vector<pair<double, ChemicalState>>& results, const vector<string>& species_names);

/**
 * @brief Utility function to compare maps of doubles.
 * 
 * @param l 
 * @param r 
 * @return true 
 * @return false 
 */
bool operator==(const map<string, double>& l, const map<string, double>& r);