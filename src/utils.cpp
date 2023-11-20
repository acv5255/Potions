#include <fstream>
#include <sstream>
#include <algorithm>
#include "potions.hpp"

using std::ofstream;
using std::stringstream;

string getOutputFilepath(string simulationName) {
    /*
        Return the filepath of the output filename for this object
        Attaches the simulation name and date-time stamp for the simulation
     */

    stringstream s;
    s << simulationName << "_results.txt";

    return s.str();
}

// bool compare_doubles(double a, double b) {
//     return  (std::abs(a - b) < 1e-12);
// }

bool SaveEquilibriumResults(const ChemicalState& chms, vector<string> species, const string& filePath) {
    try {
        ofstream file;
        file.open(filePath);
        file << "Species,Concentration\n";
        file << "name,molar\n";
        for (int i = 0; i < species.size(); i++) {
            file << species[i] << "," << chms.concentration[i] << "\n";
        }

        std::cout << "Equilibrium results written to: " << filePath  << "\n";

        file.close();

        return true;
    }
    catch (const std::exception&) {
        std::cerr << "Failed to write equilibrium results\n";
        return false;
    }
}

bool SaveKineticResults(const vector<pair<double, ChemicalState>>& res, const vector<string>& species, const string& filePath) {
    try {
        ofstream file;
        file.open(filePath);
        file << "Time,";
        for (auto x: species) file << x << ",";
        file << "\n";

        for (auto pt: res) {
            file << pt.first << ",";
            for (auto s: pt.second.concentration) file << s << ",";
            file << "\n";
        }
        file << "\n";

        file.close();
        std::cout << "Equilibrium results written to: " << filePath  << "\n";
        return true;

    }
    catch (const std::exception&) {
        std::cerr << "Failed to write kinetic results\n";
        return false;
    }
}

int getCharge(const string& name) {
    const int positiveCharge = std::count(name.cbegin(), name.cend(), '+');
    const int negativeCharge = std::count(name.cbegin(), name.cend(), '-');

    return positiveCharge - negativeCharge;
}