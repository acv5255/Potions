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
        /**
         * @brief Default constructor for ModelInputs
        */
        ModelInputs();

        /**
         * Construct the model inputs from a ChemFile object and Database object.
         * There is post-processing with some of this data to be done.
         * @brief Construct the model inputs from a ChemFile and Database
         * @param chem The ChemFile read from an 'chem.yaml' file
         * @param cdbs The Database read from the 'cdbs.yaml' file
         * @return none
        */
        ModelInputs(ChemFile chem, Database cdbs);

        /**
         * Static method for reading the input files 
         * @brief Read input files
         * @param input_name The name of the folder in the input folder containing your input files
         * @return The model inputs object after having read the inputs
        */
        static ModelInputs read_inputs(const string& inputName);
        
        /**
         * @brief Returns the the model run type, Equilibrium or Kinetic
         * @return The enum type of the model run
        */
        PotionsRunType run_type();

        /**
         * Read the initial total concentrations from the chem.yaml file and
         * initialize the total concentration vector. The concentration vector is set at a 
         * uniform value of 1e-7 moles per liter. This number is arbitrary and only intended 
         * so that the concentrations are not zero, which could be a problem later on.
         * @brief Determine initial chemical state
         * @return Returns the total concentrations at the start and a uniform concentration vector
        */
        ChemicalState initial_chem_state();

        /**
         * This function reads the list of primary and secondary species from the ChemFile object
         * and determines the stoichiometry and equilibrium parameters from the Database object. 
         * The stoichiometry matrix describes how changes in concentrations propagate throughout
         * the system, and the equilibrium parameters sets the requirement for a system at 
         * equilibrium.
         * @brief Construct the equilibrium parameters for this chemical system
         * @return The EquilibriumConstants object for this chemical system
        */
        EquilibriumConstants equilibrium_constants();
        
        /**
         * This function constructs the mass conservation matrix from the list of primary
         * species in this model. This matrix represents how each of the secondary species 
         * relate to the primary species. This model does not explicitly track H+ and instead
         * uses a charge balance for all species for the entire solution. This allows the 
         * system to have closure while also imposing another physical law on the system.
         * @brief Construct the mass conservation matrix
         * @return The TotalConstants object for this chemical system
        */
        TotalConstants total_constants();

        /**
         * This function construct the kinetic constants data structures for this chemical 
         * system. This function takes the minerals listed in the ChemFile object and 
         * reads the stoichiometry, rate constant, and kinetic equilibrium parameters from
         * the Database object. These parameters are used in the kinetic portions of the 
         * model, which uses a Transition-State-Theory-type rate law.
         * @brief Construct the kinetic mass constants
         * @return The KineticConstants object for this chemical system
        */
        KineticConstants kinetic_constants();

        /**
         * @brief Get a vector of the species names for this chemical system in the order
         * that they appear in the solution.
         * @return An ordered vector of species names 
        */
        vector<string> species_names();

        /**
         * @brief Get a map of aqueous species names and indices in the system
         * @return A map of aqueous species in this system.
        */
        map<string, unsigned int> species_map();

        /**
         * @brief Getter for the ChemFile object
         * @return The ChemFile object for this chemical system
        */
        ChemFile chem() { return _chem;};

        /**
         * @brief Getter for the Database object
         * @return The Database object for this chemical system
        */
        Database database() { return _dbs;};

        /**
         * @brief Getter for the mineral surface areas for this chemical system
         * @return The vector of mineral surface areas for this chemical system
        */
        map<string, double> surface_areas() { return _chem.mineral_surface_areas(); };
};