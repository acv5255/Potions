#include <vector>
#include <string>
#include <filesystem>
#include "model_inputs.hpp"

using std::string;
using std::filesystem::path;
using std::vector;

string getOutputFilepath(string simulationName);

bool compare_doubles(double a, double b);
bool SaveEquilibriumResults(const ChemicalState& chms, vector<string> species, const string& filePath);
bool SaveKineticResults(const vector<ChemicalState>& chms, vector<string> species, const string& filePath);
int  getCharge(const string& name);