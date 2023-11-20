#include <fstream>
#include <sstream>
#include <filesystem>
#include "potions.hpp"

using namespace std::filesystem;

bool compare_doubles(double a, double b) {
    return  (std::abs(a - b) < 1e-12);
}


PotionsRunType ModelInputs::runType() {
    return _chem.runType();
}


ModelInputs::ModelInputs() {

}


ModelInputs::ModelInputs(ChemFile chem, Database cdbs) {
    unsigned int counter = 0;
    map<string, unsigned int> speciesMap;
    for (auto x: chem.primarySpecies()) {
        speciesMap[x.first] = counter;
        counter += 1;
    }

    for (auto x: chem.secondarySpecies()) {
        speciesMap[x] = counter;
        counter += 1;
    }

    // Create the mineral map
    _mineralMap = {};
    int min_id = 0;
    for (auto x: chem.mineralSurfaceAreas()) {
        _mineralMap[x.first] = min_id;
        min_id += 1;
    }

    this->_chem = chem;
    this->_dbs = cdbs;
    this->_speciesMap = speciesMap;
}


ModelInputs ModelInputs::ReadInputs(const string& inputName) {
    // Read chem.yaml and cdbs.yaml
    const path inputDir = current_path() / "input" / inputName;
    const path cdbs_path = inputDir / "cdbs.yaml";
    const path chem_path = inputDir / "chem.yaml";

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

    ModelInputs modelInputs;

    ChemFile chem = ChemFile::FromFile(chem_path.string());
    Database cdbs = Database::FromFile(cdbs_path.string());

    // Construct the species map
    unsigned int counter = 0;
    map<string, unsigned int> speciesMap;
    for (auto x: chem.primarySpecies()) {
        speciesMap[x.first] = counter;
        counter += 1;
    }

    for (auto x: chem.secondarySpecies()) {
        speciesMap[x] = counter;
        counter += 1;
    }

    // Create the mineral map
    map<string, unsigned int> mineralMap = {};
    int min_id = 0;
    for (auto x: chem.mineralSurfaceAreas()) {
        mineralMap[x.first] = min_id;
        min_id += 1;
    }

    modelInputs._chem = chem;
    modelInputs._dbs = cdbs;
    modelInputs._speciesMap = speciesMap;
    modelInputs._mineralMap = mineralMap;

    return modelInputs;
}


ChemicalState ModelInputs::initChemState() {
    vec conc = arma::zeros(_speciesMap.size());
    conc = 1e-7;
    vec totConc = arma::zeros(_chem.primarySpecies().size());
    for (auto x: _chem.primarySpecies()) {
        if (x.first == "H+") totConc[_speciesMap.at("H+")] == 0.0;
        else totConc[_speciesMap.at(x.first)] = x.second;
    }

    return ChemicalState(conc, totConc);
}


EquilibriumConstants ModelInputs::equilibriumConstants() {
    vector<string> primSpecies;
    vector<string> secSpecies;

    for (auto x: _chem.primarySpecies()) {
        primSpecies.push_back(x.first);
    }

    for (auto x: _chem.secondarySpecies()) {
        secSpecies.push_back(x);
    }


    const int numSpecies = primSpecies.size() + secSpecies.size();
    const int numTotalSpecies = primSpecies.size();
    const int numRows = secSpecies.size();  // The number of equilibrium reactions is the number of secondary species
    const int numCols = numSpecies;

    const map<string, unsigned int> specMap = this->speciesMap();

    // Create the stoichiometry matrix
    mat stoich = arma::zeros(secSpecies.size(), numSpecies);
    vec logK = arma::zeros(secSpecies.size());

    unsigned int counter = 0;
    for (auto spec: secSpecies) {
        if (!_dbs.getSecondarySpecies().contains(spec)) {
            std::cerr << "Error: database does not contain the secondary species " << spec << "\n";
            exit(-1);
        }
        const SecondarySpecies s = _dbs.getSecondarySpecies().at(spec);

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

        logK[counter] = s.equilibriumConstant;

        counter += 1;
    }
        
    return EquilibriumConstants(stoich, logK);
}


TotalConstants ModelInputs::totalConstants() {
    const unsigned int numSpecies = _speciesMap.size();
    const unsigned int numTotalSpecies = _chem.primarySpecies().size();

    mat tot_mat = arma::zeros(numTotalSpecies, numSpecies);

    for (auto x: _chem.primarySpecies()) {
        const auto rowId = _speciesMap.at(x.first);
        const string primSpecies = x.first;

        if (x.first == "H+") {
            // Charge balance
            for (auto s: speciesNames()) {
                const auto colId = _speciesMap.at(s);
                tot_mat(rowId,colId) = getCharge(s);
            }
        }
        else {
            const auto colId = _speciesMap.at(x.first);
            // Set the primary species
            tot_mat(rowId, colId) = 1.0;
            // Go through the secondary species now
            for (auto y: _chem.secondarySpecies()) {
                if (_dbs.getSecondarySpecies().at(y).stoichiometry.contains(primSpecies)) {
                    tot_mat(rowId, _speciesMap.at(y)) = std::abs(_dbs.getSecondarySpecies().at(y).stoichiometry.at(primSpecies));
                }
            }
        }
    }

    return TotalConstants(tot_mat);
}


KineticConstants ModelInputs::kineticConstants() {
    const int numSpecies = _speciesMap.size();
    const int numMinerals = _chem.mineralSurfaceAreas().size();

    // std::cerr << "Mineral map: \n";
    // for (auto x: _mineralMap) {
    //     std::cerr << x.second << ": " << x.first << "\n";
    // }

    mat kin_mat = arma::zeros(numMinerals, numSpecies);
    vec eq_const = arma::zeros(numMinerals);
    vec kin_const = arma::zeros(numMinerals);

    for (auto chemMineral: _chem.mineralSurfaceAreas()) {
        const string mineralName = chemMineral.first;
        if (!_dbs.getMineralSpecies().contains(mineralName)) {
            std::cerr << "Database does not contain mineral: " << mineralName << "\n";
            exit(-1);
        }
        const MineralSpecies mineral = _dbs.getMineralSpecies().at(mineralName);
        const unsigned int minId = _mineralMap.at(mineralName);

        eq_const[minId] = mineral.equilibriumConstant;
        kin_const[minId] = mineral.rateConstant;

        for (auto y: mineral.stoichiometry) {
            const unsigned int colId = _speciesMap.at(y.first);
            kin_mat(minId, colId) = y.second;
        }
    }

    return KineticConstants(kin_mat, kin_const, eq_const);
}


vector<string> ModelInputs::speciesNames() {
    /*
        Return a list of the chemical species the system (in order)
     */
    vector<string> names;

    for (auto x: _chem.primarySpecies()) {
        names.push_back(x.first);
    }

    for (auto x: _chem.secondarySpecies()) {
        names.push_back(x);
    }

    return names;
}


map<string, unsigned int> ModelInputs::speciesMap() {
    return _speciesMap;
}


void ModelInputs::print() {
    // Print the model inputs
    throw NotImplemented();
}