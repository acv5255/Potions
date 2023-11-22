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
        PotionsRunType _run_type;
        map<string, double> _primary_species;
        vector<string> _secondary_species;
        map<string, double> _mineral_species;
        map<string, unsigned int> _mineral_map;
        double _end_time;
        int _num_steps;

    public:
        ChemFile();
        PotionsRunType run_type();
        static ChemFile from_file(const string& filePath);
        static ChemFile from_string(const string& yamlString);
        map<string, double> primary_species();
        vector<string> secondary_species();
        map<string, double> mineral_surface_areas();
        map<string, unsigned int> mineral_map() {
            return _mineral_map;
        }
        double endTime() {
            return _end_time;
        }
        int numSteps() {
            return _num_steps;
        }
};


// Database files
struct SecondarySpecies {
    map<string, double> stoichiometry;
    double equilibrium_constant;
    
    static SecondarySpecies from_yaml(const YAML::Node& node);
    SecondarySpecies();
    SecondarySpecies(map<string,double>& stoichiometry, double equilibriumConstant) : stoichiometry(stoichiometry), equilibrium_constant(equilibriumConstant) { };
    bool operator==(const SecondarySpecies& other);
};

bool operator==(
    const map<string, SecondarySpecies>& l,
    const map<string, SecondarySpecies>& r
);

struct MineralSpecies {
    map<string, double> stoichiometry;  // Map of the stoichiometry in this mineral
    double molar_volume;                 
    double molar_mass;          
    double rate_constant;                // Rate constant of the kinetic reaction
    double equilibrium_constant;

    static MineralSpecies from_yaml(const YAML::Node& node);
    MineralSpecies();
    MineralSpecies(
        const map<string, double>& stoichiometry,
        double molar_volume,
        double molar_mass,
        double rate_constant,
        double equilibrium_constant) :
            stoichiometry(stoichiometry),
            molar_mass(molar_mass),
            molar_volume(molar_volume),
            rate_constant(rate_constant),
            equilibrium_constant(equilibrium_constant) { };
    bool operator==(const MineralSpecies& other);
};

bool operator==(
    const map<string, MineralSpecies>& l,
    const map<string, MineralSpecies>& r);

class Database {
    private:
        vector<string> _primary_species;
        map<string, SecondarySpecies> _secondary_species;
        map<string, MineralSpecies> _mineral_kinetic;

    public:
        Database();
        Database(
            vector<string> primary_species,
            map<string, SecondarySpecies> secondary_species,
            map<string, MineralSpecies> mineral_species) :
                _primary_species(primary_species),
                _secondary_species(secondary_species),
                _mineral_kinetic(mineral_species) { }
        static Database from_file(const string& file_path);
        static Database from_string(const string& yaml_string);
        map<string, double> equilibrium_constants();
        vector<string> get_primary_species();
        map<string, SecondarySpecies> get_secondary_species();
        map<string, MineralSpecies> get_mineral_species();
};

class ModelInputs {
    private:
        map<string, unsigned int> _species_map;
        map<string, unsigned int> _mineral_map;
        PotionsRunType _run_type;
        ChemFile _chem;
        Database _dbs;
        EquilibriumConstants _equilibrium_constants;
        TotalConstants _total_constants;
        KineticConstants _kinetic_constants;

    public:
        static ModelInputs read_inputs(const string& inputName);
        ModelInputs();
        ModelInputs(ChemFile chem, Database cdbs);
        PotionsRunType run_type();
        ChemicalState initial_chem_state();
        EquilibriumConstants equilibrium_constants();
        TotalConstants total_constants();
        KineticConstants kinetic_constants();
        vector<string> species_names();
        map<string, unsigned int> species_map();
        ChemFile chem() { return _chem;};
        Database database() { return _dbs;};
        map<string, double> surface_areas() { return _chem.mineral_surface_areas(); };
};