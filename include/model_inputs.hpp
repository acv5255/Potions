#pragma once

#include <filesystem>
#include <string>
#include <map>
#include <vector>
#include "yaml-cpp/yaml.h"

using std::filesystem::path;
using std::string;
using std::map;
using std::vector;

class ChemFile {
    map<string, double> _primarySpecies;
    vector<string> _secondarySpecies;

    public:
        static ChemFile FromFile(const path& filePath);
        map<string, double> primarySpecies();
        vector<string> secondarySpecies();
};

struct SecondaryEquilibriumEntry {
    string name;
    map<string, double> stoichiometry;
    double equilibriumConstant;
    
    static SecondaryEquilibriumEntry FromYaml(const YAML::Node& node);
};

struct MineralKineticEntry {
    string name;                        // Name of the mineral
    map<string, double> stoichiometry;  // Map of the stoichiometry in this mineral
    double molarVolume;                 
    double molarMass;          
    double rateConstant;                // Rate constant of the kinetic reaction
    double equilibriumConstant;

    static MineralKineticEntry FromYaml(const YAML::Node& node);
};

class Database {
    vector<string> primarySpecies;
    map<string, double> equilibriumConstants;


    public:
        static Database FromFile(const path& filePath);
        map<string, double> getEquilibriumConstants();
        vector<string> getPrimarySpecies();
};

class ModelInputs {
    private:
        PotionsRunType _runType;

    public:
        static ModelInputs ReadInputs(const string& inputName);
        PotionsRunType runType();
};