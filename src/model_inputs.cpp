#include <fstream>
#include <sstream>
#include <filesystem>
#include <thread>
#include <future>
#include "potions.hpp"

using namespace std::filesystem;

/**
 * This function is used to compare floating point numbers in a safe manner.
 * I set the error tolerance at 1e-12 becuase it is smaller than any number 
 * I should need to encounter in this problem, but this could easily be changed.
 * @brief Compare floating point numbers
 * @param a first value
 * @param b second value
 * @return Whether the values are equal to a tolerance of 1e-12
*/
bool compare_doubles(double a, double b) {
    return  (std::abs(a - b) < 1e-12);
}


PotionsRunType ModelInputs::run_type() {
    return _chem.run_type();
}



ModelInputs::ModelInputs() {

}



ModelInputs::ModelInputs(ChemFile chem, Database cdbs) {
    unsigned int counter = 0;
    map<string, unsigned int> speciesMap;
    for (auto x: chem.primary_species()) {
        speciesMap[x.first] = counter;
        counter += 1;
    }

    for (auto x: chem.secondary_species()) {
        speciesMap[x] = counter;
        counter += 1;
    }

    // Create the mineral map
    _mineral_map = {};
    int min_id = 0;
    for (auto x: chem.mineral_surface_areas()) {
        _mineral_map[x.first] = min_id;
        min_id += 1;
    }

    this->_chem = chem;
    this->_dbs = cdbs;
    this->_species_map = speciesMap;
}


ModelInputs ModelInputs::read_inputs(const string& input_name) {
    // Read chem.yaml and cdbs.yaml
    const path inputDir = current_path() / "input" / input_name;
    const path cdbs_path = inputDir / "cdbs.yaml";
    const path chem_path = inputDir / "chem.yaml";

    // Check to make sure all required directories exist
    if (!exists(inputDir)) {
        std::cerr << "Input directory at " << inputDir << " does not exist\n";
        exit(-1);
    }

    if (!exists(cdbs_path)) {
        std::cerr << "Database at " << cdbs_path << " does not exist\n";
        exit(-1);
    }
    
    if (!exists(chem_path)) {
        std::cerr << "Chemistry file at " << chem_path << " does not exist\n";
        exit(-1);
    }

    ModelInputs model_inputs;

    // For large database and chem files, construct these data structure asynchronously
    std::future<ChemFile> chem_p = std::async(ChemFile::from_file, chem_path.string());
    ChemFile chem = chem_p.get();

    std::future<Database> cdbs_p = std::async(Database::from_file, cdbs_path.string());
    Database cdbs = cdbs_p.get();

    // Construct the species map
    unsigned int counter = 0;
    map<string, unsigned int> species_map;
    for (auto x: chem.primary_species()) {
        species_map[x.first] = counter;
        counter += 1;
    }

    for (auto x: chem.secondary_species()) {
        species_map[x] = counter;
        counter += 1;
    }

    // Create the mineral map
    map<string, unsigned int> mineral_map = {};
    int min_id = 0;
    for (auto x: chem.mineral_surface_areas()) {
        mineral_map[x.first] = min_id;
        min_id += 1;
    }

    model_inputs._chem = chem;
    model_inputs._dbs = cdbs;
    model_inputs._species_map = species_map;
    model_inputs._mineral_map = mineral_map;

    return model_inputs;
}


ChemicalState ModelInputs::initial_chem_state() {
    vec conc = arma::zeros(_species_map.size());
    conc = 1e-7;
    vec totConc = arma::zeros(_chem.primary_species().size());
    for (auto x: _chem.primary_species()) {
        if (x.first == "H+") totConc[_species_map.at("H+")] == 0.0;
        else totConc[_species_map.at(x.first)] = x.second;
    }

    return ChemicalState(conc, totConc);
}


EquilibriumConstants ModelInputs::equilibrium_constants() {
    vector<string> primSpecies;
    vector<string> secSpecies;

    for (auto x: _chem.primary_species()) {
        primSpecies.push_back(x.first);
    }

    for (auto x: _chem.secondary_species()) {
        secSpecies.push_back(x);
    }


    const int numSpecies = primSpecies.size() + secSpecies.size();
    const int numTotalSpecies = primSpecies.size();
    const int numRows = secSpecies.size();  // The number of equilibrium reactions is the number of secondary species
    const int numCols = numSpecies;

    const map<string, unsigned int> specMap = this->species_map();

    // Create the stoichiometry matrix
    mat stoich = arma::zeros(secSpecies.size(), numSpecies);
    vec logK = arma::zeros(secSpecies.size());

    unsigned int counter = 0;
    for (auto spec: secSpecies) {
        if (!_dbs.get_secondary_species().contains(spec)) {
            std::cerr << "Error: database does not contain the secondary species " << spec << "\n";
            exit(-1);
        }
        const SecondarySpecies s = _dbs.get_secondary_species().at(spec);

        if (!specMap.contains(spec)) {
            std::cerr << "Species map does not contain entry: '" << spec << "'\n";
            exit(-1);
        }
        stoich(counter, specMap.at(spec)) = 1.0;    // Set the value for this species

        // Set the secondary species values
        for (auto entry: s.stoichiometry) {
            const unsigned int index = specMap.at( entry.first );
            stoich(counter, index ) = entry.second;
        }

        logK[counter] = s.log_eq_const;

        counter += 1;
    }
        
    return EquilibriumConstants(stoich, logK);
}


TotalConstants ModelInputs::total_constants() {
    const unsigned int numSpecies = _species_map.size();
    const unsigned int numTotalSpecies = _chem.primary_species().size();

    mat tot_mat = arma::zeros(numTotalSpecies, numSpecies);

    for (auto x: _chem.primary_species()) {
        const auto rowId = _species_map.at(x.first);
        const string primSpecies = x.first;

        if (x.first == "H+") {
            // Charge balance
            for (auto s: species_names()) {
                const auto colId = _species_map.at(s);
                tot_mat(rowId,colId) = get_charge(s);
            }
        }
        else {
            const auto colId = _species_map.at(x.first);
            // Set the primary species
            tot_mat(rowId, colId) = 1.0;
            // Go through the secondary species now
            for (auto y: _chem.secondary_species()) {
                if (_dbs.get_secondary_species().at(y).stoichiometry.contains(primSpecies)) {
                    tot_mat(rowId, _species_map.at(y)) = std::abs(_dbs.get_secondary_species().at(y).stoichiometry.at(primSpecies));
                }
            }
        }
    }

    return TotalConstants(tot_mat);
}


KineticConstants ModelInputs::kinetic_constants() {
    const int numSpecies = _species_map.size();
    const int numMinerals = _chem.mineral_surface_areas().size();

    mat kin_mat = arma::zeros(numMinerals, numSpecies);
    vec eq_const = arma::zeros(numMinerals);
    vec kin_const = arma::zeros(numMinerals);

    for (auto chemMineral: _chem.mineral_surface_areas()) {
        const string mineralName = chemMineral.first;
        if (!_dbs.get_mineral_species().contains(mineralName)) {
            std::cerr << "Database does not contain mineral: " << mineralName << "\n";
            exit(-1);
        }
        const MineralSpecies mineral = _dbs.get_mineral_species().at(mineralName);
        const unsigned int minId = _mineral_map.at(mineralName);

        eq_const[minId] = mineral.log_eq_const;
        kin_const[minId] = mineral.log_kin_const;

        for (auto y: mineral.stoichiometry) {
            const unsigned int colId = _species_map.at(y.first);
            kin_mat(minId, colId) = y.second;
        }
    }

    return KineticConstants(kin_mat, kin_const, eq_const);
}


vector<string> ModelInputs::species_names() {
    /*
        Return a list of the chemical species the system (in order)
     */
    vector<string> names;

    for (auto x: _chem.primary_species()) {
        names.push_back(x.first);
    }

    for (auto x: _chem.secondary_species()) {
        names.push_back(x);
    }

    return names;
}


map<string, unsigned int> ModelInputs::species_map() {
    return _species_map;
}
