#pragma once

/** @file */

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
const string SIMULATION_TYPE_TOKEN = "SimulationType";  /** Name of the parameter for species */
const string PRIMARY_SECTION_TOKEN = "PrimarySpecies"; /** Token for the primary species block in 'chem.yaml' */
const string SECONDARY_SECTION_TOKEN = "SecondarySpecies"; /** Token for secondary species block in 'chem.yaml' */
const string MINERAL_SECTION_TOKEN = "MineralSpecies"; /** Token for mineral species block in 'chem.yaml' */
const string SIMTYPE_EQ_TOKEN = "Equilibrium";  /** Token specifying an equilibrium simulation in 'chem.yaml' */
const string SIMTYPE_KIN_TOKEN = "Kinetic";     /** Token specifying a kinetic simulation in 'chem.yaml' */
const string END_TIME_TOKEN = "EndTime";        /** Token specifying the time to run the simulation to */
const string NUM_STEPS = "NumSteps";            /** Token for the number of steps to use in the simulation */

// Database tokens type
const string DBS_ID_TOKEN = "Database";         /** Token for the database ID token */
const string DBS_PRIMARY_SECTION_TOKEN = "PrimarySpecies";  /** Token specifying the primary section list in 'cdbs.yaml' */
const string DBS_SECONDARY_SECTION_TOKEN = "SecondarySpecies"; /** Token specifying the secondary block in 'cdbs.yaml' */
const string DBS_MINERAL_SECTION_TOKEN = "MineralSpecies";  /** Token specifying the mineral section in 'cdbs.yaml' */
const string DBS_STOICHIOMETRY_TOKEN = "stoichiometry";     /** Token specifying a stoichiometry block in the database */
const string DBS_LOGK_TOKEN = "logK";                       /** Token specifying the equilibrium constant in the database */
const string DBS_MOLAR_MASS_TOKEN = "molarMass";            /** Token specifying the molar mass the 'cdbs.yaml' */
const string DBS_MOLAR_VOLUME_TOKEN = "molarVolume";        /** Token specifying the molar volume in the 'cdbs.yaml' */
const string DBS_RATE_CONST_TOKEN = "rate";                 /** Token specifying the kinetic rate constant in 'cdbs.yaml' */


/**
 * @brief Chemical system parameters
 * The class that contains the complete information about the primary species.
 * This file contains the list of primary species including the total concentrations
 * for each primary species. The PrimarySpecies block requires both the name of the
 * primary species and the total concentration. This section also requires the
 * secondary species block to be populated with at least 1 secondary species. For a kinetic
 * simulation, this requires the MineralKinetic block to have at least one mineral entry
 * which contains both a mineral name and a mineral surface area.
 */
class ChemFile {

    private:
        PotionsRunType _run_type;               /** Whether the run is an equilibrium or kinetic solution */
        map<string, double> _primary_species;   /** Map of primary species total concentrations */
        vector<string> _secondary_species;      /** List of secondary species in the system*/
        map<string, double> _mineral_species;   /** Map of mineral surface areas in the*/
        map<string, unsigned int> _mineral_map; /** Map of the mineral indices */
        double _end_time;                       /** Time (in seconds) to end the simulation */
        int _num_steps;                         /** Number of steps to run the simulation */

    public:
        ChemFile();                             /** Default constructor */
        PotionsRunType run_type();              /** Enum value for whether this is an equilibrium or kinetic simulation*/

        /**
         * @brief Load 'chem.yaml' from a file
         * Read a 'chem.yaml' file and produce a ChemFile object. This function exits if the 
         * 'chem.yaml' file is not formatted correctly.
         * @param file_path The file path to the 'chem.yaml' file.
         */
        static ChemFile from_file(const string& file_path);

        /**
         * @brief Produce a ChemFile object from a string containing the contents of a 'chem.yaml' file.
         * @param yaml_string A string containing the entire contents of 'chem.yaml'
         * @return ChemFile
         */
        static ChemFile from_string(const string& yaml_string);

        /**
         * @brief Primary species map getter
         * @return map<string, double> 
         */
        map<string, double> primary_species();

        /**
         * @brief Secondary species vector getter
         * 
         * @return vector<string> 
         */
        vector<string> secondary_species();
        
        /**
         * @brief Mineral surface area vector getter
         * 
         * @return map<string, double> 
         */
        map<string, double> mineral_surface_areas();

        /**
         * @brief Mineral index map getter
         * 
         * @return map<string, unsigned int> 
         */
        map<string, unsigned int> mineral_map() {
            return _mineral_map;
        }

        /**
         * @brief End time getter
         * 
         * @return double The final time (in seconds) of a kinetic simulation
         */
        double end_time() {
            return _end_time;
        }

        /**
         * @brief Getter for the number of steps in the kinetic simulation
         * 
         * @return int The number of steps for a kinetic simulation
         */
        int num_steps() {
            return _num_steps;
        }
};


/**
 * @brief The secondary species information struct
 * This data structure contains the information for how to represent a secondary species
 * in terms of primary species.
 */
struct SecondarySpecies {
    map<string, double> stoichiometry;  /** Map of stoichiometric coefficients of primary species */
    double log_eq_const;        /** The base-10 logarithm of the equilibrium constant */
    
    /**
     * @brief Produce a SecondarySpecies object from a YAML node object
     * 
     * @param node The YAML node object representing the secondary species object
     * @return SecondarySpecies 
     */
    static SecondarySpecies from_yaml(const YAML::Node& node);
    
    /**
     * @brief Construct a new Secondary Species object
     * 
     */
    SecondarySpecies();

    /**
     * @brief Construct a new Secondary Species object
     * 
     * @param stoichiometry The stoichiometric coefficients for this object
     * @param log_eq_const The base-10 logarithm of the equilibrium constant for this species
     */
    SecondarySpecies(map<string,double>& stoichiometry, double log_eq_const) : stoichiometry(stoichiometry), log_eq_const(log_eq_const) { };
    
    /**
     * @brief Compare two SecondarySpecies objects
     * 
     * @param other 
     * @return true 
     * @return false 
     */
    bool operator==(const SecondarySpecies& other);
};

/**
 * @brief Compare two maps of secondary species objects
 * 
 * @param l 
 * @param r 
 * @return true 
 * @return false 
 */
bool operator==(
    const map<string, SecondarySpecies>& l,
    const map<string, SecondarySpecies>& r
);

/**
 * @brief Data structure containing information for a mineral species
 * 
 */
struct MineralSpecies {
    map<string, double> stoichiometry;  /** Map for stoichiometry of the dissolution of this mineral */
    double molar_volume;                /** Volume per mole of mineral, in cm^3/mol */ 
    double molar_mass;                  /** Mass of a mole of mineral, in g/mol */
    double log_kin_const;               /** Base-10 logarithm of the kinetic constant */
    double log_eq_const;                /** Base-10 logarithm of the equilibrium constant */

    /**
     * @brief Produce a MineralSpecies object from a YAML node object
     * 
     * @param node 
     * @return MineralSpecies 
     */
    static MineralSpecies from_yaml(const YAML::Node& node);

    /**
     * @brief Construct a new Mineral Species object
     * 
     */
    MineralSpecies();

    /**
     * @brief Construct a new Mineral Species object
     * 
     * @param stoichiometry 
     * @param molar_volume 
     * @param molar_mass 
     * @param rate_constant 
     * @param equilibrium_constant 
     */
    MineralSpecies(
        const map<string, double>& stoichiometry,
        double molar_volume,
        double molar_mass,
        double rate_constant,
        double equilibrium_constant) :
            stoichiometry(stoichiometry),
            molar_mass(molar_mass),
            molar_volume(molar_volume),
            log_kin_const(rate_constant),
            log_eq_const(equilibrium_constant) { }

    /**
     * @brief Compare two MineralSpecies objects
     * 
     * @param other 
     * @return true 
     * @return false 
     */
    bool operator==(const MineralSpecies& other);
};

/**
 * @brief Compare maps of MineralSpecies objects
 * 
 * @param l 
 * @param r 
 * @return true 
 * @return false 
 */
bool operator==(
    const map<string, MineralSpecies>& l,
    const map<string, MineralSpecies>& r);


/**
 * @brief Data structure for data contained in the chemical database, 'cdbs.yaml'
 * This data structure contains all of the chemical information that the user does
 * not need to change. In the database, there is a list of primary species, 
 * secondary species, including the equilibrium parameters, and mineral species 
 * including the kinetic and equilibrium parameters.
 */
class Database {
    private:
        vector<string> _primary_species;                    /** List of primary species in the database */
        map<string, SecondarySpecies> _secondary_species;   /** Map of secondary species information in the database */
        map<string, MineralSpecies> _mineral_kinetic;       /** Map of mineral species information in the database */

    public:
        /**
         * @brief Construct a new Database object
         * 
         */
        Database();

        /**
         * @brief Construct a new Database object
         * 
         * @param primary_species 
         * @param secondary_species 
         * @param mineral_species 
         */
        Database(
            vector<string> primary_species,
            map<string, SecondarySpecies> secondary_species,
            map<string, MineralSpecies> mineral_species) :
                _primary_species(primary_species),
                _secondary_species(secondary_species),
                _mineral_kinetic(mineral_species) { }

        /**
         * @brief Read a 'cdbs.yaml' file
         * 
         * @param file_path The file path to the 'cdbs.yaml'
         * @return Database 
         */
        static Database from_file(const string& file_path);

        /**
         * @brief Produce a Database object from a string containing the contents of the 'cdbs.yaml' file
         * 
         * @param yaml_string 
         * @return Database 
         */
        static Database from_string(const string& yaml_string);

        /**
         * @brief Get the primary species vector
         * 
         * @return vector<string> 
         */
        vector<string> get_primary_species();

        /**
         * @brief Get the secondary species map
         * 
         * @return map<string, SecondarySpecies> 
         */
        map<string, SecondarySpecies> get_secondary_species();

        /**
         * @brief Get the mineral kinetic species map
         * 
         * @return map<string, MineralSpecies> 
         */
        map<string, MineralSpecies> get_mineral_species();
};

/**
 * @brief Wrapper for the model inputs
 * Data structure containing the input data to run the simulations.
 * This is a way to interface with both the ChemFile and Database objects
 * to simplify the model interface, as this class contains methods to
 * construct the mass conservation, equilibrium, and kinetic constants.
 */
class ModelInputs {
    private:
        map<string, unsigned int> _species_map;         /** Indices of the chemical species in the concentration vector */
        map<string, unsigned int> _mineral_map;         /** Indices of the minerals in the kinetic map */
        PotionsRunType _run_type;                       /** The model run type - Kinetic or Equilibrium */
        ChemFile _chem;                                 /** The ChemFile object for the simulation */
        Database _dbs;                                  /** The Database object for the simulation */

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