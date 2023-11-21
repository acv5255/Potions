#include <fstream>
#include <sstream>
#include <algorithm>
#include "potions.hpp"
#include <matplot/matplot.h>

using std::ofstream;
using std::stringstream;
using namespace matplot;

string getOutputFilepath(string simulationName) {
    /*
        Return the filepath of the output filename for this object
        Attaches the simulation name and date-time stamp for the simulation
     */

    stringstream s;
    s << simulationName << "_results.txt";

    return s.str();
}


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


template<typename T>
void PrintMatrix(const Mat<T>& m) {
    for (int i = 0; i < m.n_rows; i++) {
        for (int j = 0; j < m.n_cols; j++) {
            std::cout << m(i,j) << " ";
        }

        std::cout << "\n";
    }
}

bool PlotResults(const vector<pair<double, ChemicalState>>& results, const vector<string>& speciesNames) {
    /*
        Plot the kinetic time-series results
     */
    // Construct input arrays 
    const int N = results.size();
    const int numSpecies = speciesNames.size();
    vector<double> ts(N);
    vector<vector<double>> data(numSpecies);
    for (int i = 0; i < numSpecies; i++) {
        data[i] = vector<double>(N);
    }

    int counter = 0;
    for (auto x: results) {
        ts[counter] = x.first;
        for (int j = 0; j < numSpecies; j++) {
            data[j][counter] = x.second.concentration[j];
        }

        counter += 1;
    }

    // Plot results
    hold(on);
    for (int i = 0; i < numSpecies; i++) {
        semilogy(ts, data[i]);
    }
    legend(speciesNames);

    // Show plot
    show();
    title("Kinetic solution");
    xlabel("Time [s]");
    ylabel("Concentration [molar]");

    return true;
}