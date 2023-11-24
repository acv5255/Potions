#include "potions.hpp"

// Default constructor, never used
ChemFile::ChemFile() {

}

/*
    Read the chem.yaml input file and produce a chemfile object.
 */
ChemFile ChemFile::from_file(const string& file_path) {
    // Read the file from a string into a ChemFile object
    std::ifstream file(file_path);
    std::stringstream buffer;
    buffer << file.rdbuf();

    return ChemFile::from_string(buffer.str());
}

/*
    Produce a ChemFile object from a complete file string.
    This function exists to make testing easier.
 */
ChemFile ChemFile::from_string(const string& yaml_string) {
    YAML::Node root = YAML::Load(yaml_string);

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

    // End time
    double endTime = 0.0;
    if (root[END_TIME_TOKEN]) {
        endTime = root[END_TIME_TOKEN].as<double>();

        if (endTime < 0.0) {
            std::cerr << "End time for kinetic simulation must be > 0\n";
            exit(-1);
        }
    }

    // Number of simulation steps
    int numSteps = 1;
    if (root[NUM_STEPS]) {
        numSteps = root[NUM_STEPS].as<int>();

        if (numSteps < 1) {
            std::cerr << "Number of steps must be > 0\n";
            exit(-1);
        }
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

    map<string, double> minAreas;       // Mineral surface areas
    if (minNode.size() > 0) {
        if (!minNode.IsMap()) {
            std::cerr << "Error: Mineral species section must be a YAML Map\n";
            exit(-1);
        }
        for (YAML::const_iterator it=minNode.begin(); it != minNode.end(); it++) {
            minAreas[it->first.as<string>()] = it->second.as<double>();
        }
    }

    map<string, unsigned int> minMap;
    unsigned int count = 0;
    for (auto x: minAreas) {
        minMap[x.first] = count;
        count += 1;
    }
    //=====================================//

    // Construct object and return
    ChemFile chem;
    chem._run_type = simulationType;
    chem._primary_species = primSpec;
    chem._secondary_species = secSpec;
    chem._mineral_species = minAreas;
    chem._end_time = endTime;
    chem._num_steps = numSteps;
    chem._mineral_map = minMap;

    return chem;
}

map<string, double> ChemFile::primary_species() {
    return _primary_species;
}

vector<string> ChemFile::secondary_species() {
    return _secondary_species;
}

map<string, double> ChemFile::mineral_surface_areas() {
    return _mineral_species;
}

PotionsRunType ChemFile::run_type() {
    return _run_type;
}