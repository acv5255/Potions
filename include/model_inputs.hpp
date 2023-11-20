#pragma once

#include <filesystem>
#include <string>
#include <map>
#include <vector>
#include "yaml-cpp/yaml.h"
#include "potions.hpp"

using std::filesystem::path;
using std::string;
using std::map;
using std::vector;

// Constants //
// Chem file tokens
const string SIMULATION_TYPE_TOKEN = "SimulationType";
const string PRIMARY_SECTION_TOKEN = "PrimarySpecies";
const string SECONDARY_SECTION_TOKEN = "SecondarySpecies";
const string MINERAL_SECTION_TOKEN = "MineralSpecies";
const string SIMTYPE_EQ_TOKEN = "Equilibrium";  // Specifies an equilibrium simulation
const string SIMTYPE_KIN_TOKEN = "Kinetic";     // Specifies a kinetic simulation
const string END_TIME_TOKEN = "EndTime";
const string NUM_STEPS = "NumSteps";

// Database tokens type
const string DBS_ID_TOKEN = "Database";
const string DBS_PRIMARY_SECTION_TOKEN = "PrimarySpecies";
const string DBS_SECONDARY_SECTION_TOKEN = "SecondarySpecies";
const string DBS_MINERAL_SECTION_TOKEN = "MineralSpecies";
const string DBS_STOICHIOMETRY_TOKEN = "stoichiometry";
const string DBS_LOGK_TOKEN = "logK";
const string DBS_MOLAR_MASS_TOKEN = "molarMass";
const string DBS_MOLAR_VOLUME_TOKEN = "molarVolume";
const string DBS_RATE_CONST_TOKEN = "rate";

bool operator==(const map<string, double>& l, const map<string, double>& r);

class ChemFile {
    private:
        PotionsRunType _runType;
        map<string, double> _primarySpecies;
        vector<string> _secondarySpecies;
        map<string, double> _mineralSpecies;
        map<string, unsigned int> _mineralMap;
        double _endTime;
        int _numSteps;

    public:
        ChemFile();
        PotionsRunType runType();
        static ChemFile FromFile(const string& filePath);
        static ChemFile FromString(const string& yamlString);
        map<string, double> primarySpecies();
        vector<string> secondarySpecies();
        map<string, double> mineralSurfaceAreas();
        map<string, unsigned int> mineralMap() { return _mineralMap; };
        double endTime() { return _endTime; };
        int numSteps() { return _numSteps; };
};

// Database files
struct SecondarySpecies {
    map<string, double> stoichiometry;
    double equilibriumConstant;
    
    static SecondarySpecies FromYaml(const YAML::Node& node);
    SecondarySpecies();
    SecondarySpecies(map<string,double>& stoichiometry, double equilibriumConstant) : stoichiometry(stoichiometry), equilibriumConstant(equilibriumConstant) { };
    bool operator==(const SecondarySpecies& other);
};

bool operator==(const map<string, SecondarySpecies>& l, const map<string, SecondarySpecies>& r);

struct MineralSpecies {
    map<string, double> stoichiometry;  // Map of the stoichiometry in this mineral
    double molarVolume;                 
    double molarMass;          
    double rateConstant;                // Rate constant of the kinetic reaction
    double equilibriumConstant;

    static MineralSpecies FromYaml(const YAML::Node& node);
    MineralSpecies();
    MineralSpecies(
        const map<string, double>& stoichiometry, 
        double molarVolume, 
        double molarMass,
        double rateConstant,
        double equilibriumConstant
        ) : stoichiometry(stoichiometry), molarMass(molarMass), molarVolume(molarVolume),
            rateConstant(rateConstant), equilibriumConstant(equilibriumConstant) { };
    bool operator==(const MineralSpecies& other);
};

bool operator==(const map<string, MineralSpecies>& l, const map<string, MineralSpecies>& r);

class Database {
    private:
        vector<string> _primSpec;
        map<string, SecondarySpecies> _secSpec;
        map<string, MineralSpecies> _minKin;

    public:
        Database();
        Database(vector<string> primSpecies, map<string, SecondarySpecies> secSpecies, map<string, MineralSpecies> minSpecies) :
            _primSpec(primSpecies), _secSpec(secSpecies), _minKin(minSpecies) { };
        static Database FromFile(const string& filePath);
        static Database FromString(const string& yamlString);
        map<string, double> getEquilibriumConstants();
        vector<string> getPrimarySpecies();
        map<string, SecondarySpecies> getSecondarySpecies();
        map<string, MineralSpecies> getMineralSpecies();
};

class ModelInputs {
    private:
        map<string, unsigned int> _speciesMap;
        map<string, unsigned int> _mineralMap;
        PotionsRunType _runType;
        ChemFile _chem;
        Database _dbs;
        EquilibriumConstants _eq_consts;
        TotalConstants _tot_consts;
        KineticConstants _kin_consts;

    public:
        static ModelInputs ReadInputs(const string& inputName);
        ModelInputs();
        ModelInputs(ChemFile chem, Database cdbs);
        PotionsRunType runType();
        ChemicalState initChemState();
        EquilibriumConstants equilibriumConstants();
        TotalConstants totalConstants();
        KineticConstants kineticConstants();
        vector<string> speciesNames();
        map<string, unsigned int> speciesMap();
        ChemFile chem() { return _chem;};
        Database database() { return _dbs;};
        map<string, double> surfaceAreas() { return _chem.mineralSurfaceAreas(); };
        // map<string, unsigned int> mineralMap() { return _chem.min}
        void print();
};