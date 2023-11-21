#include <vector>
#include <string>
#include <filesystem>
#include <utility>
#include "model_inputs.hpp"

using std::string;
using std::filesystem::path;
using std::vector;
using std::pair;

string getOutputFilepath(string simulationName);

bool compare_doubles(double a, double b);
bool SaveEquilibriumResults(const ChemicalState& chms, vector<string> species, const string& filePath);
bool SaveKineticResults(const vector<pair<double, ChemicalState>>& res, const vector<string>& species, const string& filePath);
int  getCharge(const string& name);

template<typename T>
void PrintMatrix(const Mat<T>& m);

bool PlotResults(const vector<pair<double, ChemicalState>>& results, const vector<string>& speciesNames);