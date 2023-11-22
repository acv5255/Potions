#include <vector>
#include <string>
#include <filesystem>
#include <utility>
#include "model_inputs.hpp"

using std::string;
using std::filesystem::path;
using std::vector;
using std::pair;

string get_output_file_path(string simulation_name);

bool compare_doubles(double a, double b);
bool save_equilibrium_results(const ChemicalState& chms, vector<string> species, const string& file_path);
bool save_kinetic_results(const vector<pair<double, ChemicalState>>& res, const vector<string>& species, const string& filePath);
int  get_charge(const string& name);

template<typename T>
void print_matrix(const Mat<T>& m);

bool plot_results(const vector<pair<double, ChemicalState>>& results, const vector<string>& speciesNames);