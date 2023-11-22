#include <fstream>
#include <sstream>
#include <algorithm>
#include <chrono>
#include <ctime>
#include <string.h>
#include <array>
#include "potions.hpp"
#include <matplot/matplot.h>

using std::ofstream;
using std::stringstream;
namespace plt = matplot;
using namespace std::chrono;
using std::time_t;
using std::ctime;
using std::array;

/*
    Return a timestamp in the format of YYYYMMHH_HHMMSS
 */
string get_timestamp() {
    const time_point now = system_clock::now();
    const time_t end_time = system_clock::to_time_t(now);
    tm *cur_time = std::localtime(&end_time);

    array<char, 1024> ts;

    strftime(ts.data(), ts.size(), "%Y%m%d_%H%M%S", cur_time);

    return string(ts.data());
}


/*
    Construct a simple output file in the format:
    'simname_timestamp_results.csv'
*/
string get_output_file_path(string simulation_name) {
    /*
        Return the filepath of the output filename for this object
        Attaches the simulation name and date-time stamp for the simulation
     */

    stringstream s;
    s << simulation_name << "_" << get_timestamp() << "_results.csv";

    return s.str();
}


/*
    Write the results of an equilibrium simulation to an output file
 */
bool save_equilibrium_results(const ChemicalState& chms, vector<string> species, const string& filePath) {
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


/*
    Write the kinetic simulation to an CSV file
 */
bool save_kinetic_results(const vector<pair<double, ChemicalState>>& res, const vector<string>& species, const string& filePath) {
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


/*
    Determine the ionic charge of a chemical species based on the name.
    Example:
        'Ca++', calcium ion, has an ionic charge of 2
        'Cl-', chloride ion, has an ionic charge of -1
        'H2CO3', carbonic acid, has an ionic charge of 0
 */
int get_charge(const string& name) {
    const int positiveCharge = std::count(name.cbegin(), name.cend(), '+');
    const int negativeCharge = std::count(name.cbegin(), name.cend(), '-');

    return positiveCharge - negativeCharge;
}


template<typename T>
void print_matrix(const Mat<T>& m) {
    for (int i = 0; i < m.n_rows; i++) {
        for (int j = 0; j < m.n_cols; j++) {
            std::cout << m(i,j) << " ";
        }

        std::cout << "\n";
    }
}


/*
    Plot the kinetic solution in a GUI with plots that can be saved
    to a file.
*/
bool plot_results(const vector<pair<double, ChemicalState>>& results, const vector<string>& speciesNames) {
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
    plt::hold(plt::on);
    for (int i = 0; i < numSpecies; i++) {
        plt::semilogy(ts, data[i]);
    }
    plt::legend(speciesNames);

    // Show plot
    plt::title("Kinetic solution");
    plt::xlabel("Time [s]");
    plt::ylabel("Concentration [molar]");
    plt::show();

    return true;
}