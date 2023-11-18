#include "potions.hpp"

ChemFile::ChemFile() {

}

ChemFile ChemFile::FromFile(const string& filePath) {
    // Read the file from a string into a ChemFile object
    std::ifstream file(filePath);
    std::stringstream buffer;
    buffer << file.rdbuf();

    return ChemFile::FromString(buffer.str());
}

ChemFile ChemFile::FromString(const string& yamlString) {
    YAML::Node root = YAML::Load(yamlString);

    // Check for required sections
    if (!root[SIMULATION_TYPE_TOKEN]) {
        std::cerr << "Error: Missing primary species section in 'chem.yaml'\n";
        exit(-1);
    }

    if (!root[PRIMARY_SECTION_TOKEN]) {
        std::cerr << "Error: Missing primary species section in 'chem.yaml'\n";
        exit(-1);
    }

    if (!root[SECONDARY_SECTION_TOKEN]) {
        std::cerr << "Error: Missing primary species section in 'chem.yaml'\n";
        exit(-1);
    }

    if (!root[MINERAL_SECTION_TOKEN]) {
        std::cerr << "Error: Missing primary species section in 'chem.yaml'\n";
        exit(-1);
    }

    const YAML::Node primNode = root[PRIMARY_SECTION_TOKEN];
    const YAML::Node secNode = root[SECONDARY_SECTION_TOKEN];
    const YAML::Node minNode = root[MINERAL_SECTION_TOKEN];

    //===== Now, read in each section =====//
    // Simulation type
    PotionsRunType simulationType;
    const string simTypeName = root[SIMULATION_TYPE_TOKEN].as<string>();
    if (simTypeName == SIMTYPE_KIN_TOKEN) {
        simulationType = PotionsRunType::KINETIC;
    }
    else if (simTypeName == SIMTYPE_EQ_TOKEN) {
        simulationType = PotionsRunType::EQUILIBRIUM;
    }
    else {
        std::cerr << "Unknown simulation type: " << simTypeName << "\n";
        exit(-1);
    }

    // Primary species
    if (!(primNode.Type() == YAML::NodeType::Map)) {
        std::cerr << "Error: Primary species section must be a YAML Map\n";
        exit(-1);
    }

    map<string, double> primSpec;   // Primary species
    for (YAML::const_iterator it=primNode.begin(); it != primNode.end(); it++) {
        primSpec[it->first.as<string>()] = it->second.as<double>();
    }

    // Secondary species
    if (!secNode.IsSequence()) {
        std::cerr << "Error: Secondary species section must be a YAML Sequence\n";
        exit(-1);
    }

    vector<string> secSpec;
    for (YAML::const_iterator it=secNode.begin(); it != secNode.end(); it++) {
        secSpec.push_back(it->as<string>());
    }

    // Mineral species
    if (!minNode.IsMap()) {
        std::cerr << "Error: Mineral species section must be a YAML Map\n";
        exit(-1);
    }

    map<string, double> minAreas;       // Mineral surface areas
    for (YAML::const_iterator it=minNode.begin(); it != minNode.end(); it++) {
        minAreas[it->first.as<string>()] = it->second.as<double>();
    }
    //=====================================//

    // Construct object and return
    ChemFile chem;
    chem._runType = simulationType;
    chem._primarySpecies = primSpec;
    chem._secondarySpecies = secSpec;
    chem._mineralSpecies = minAreas;

    return chem;
}

map<string, double> ChemFile::primarySpecies() {
    return _primarySpecies;
}

vector<string> ChemFile::secondarySpecies() {
    return _secondarySpecies;
}

map<string, double> ChemFile::mineralSpecies() {
    return _mineralSpecies;
}

PotionsRunType ChemFile::runType() {
    return _runType;
}