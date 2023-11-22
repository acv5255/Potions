#include "potions.hpp"

Database::Database() {

}

Database Database::from_file(const string& filePath) {
    std::ifstream file(filePath);
    std::stringstream buffer;
    buffer << file.rdbuf();

    return Database::from_string(buffer.str());
}

bool operator==(const map<string, double>& a, const map<string, double>& b) {
    if (a.size() != b.size()) return false;

    for (auto a_i: a) {
        if (!b.contains(a_i.first)) {
            return false;
        }

        if (!(std::abs(a_i.second - b.at(a_i.first)) - 1e-12)) {
        return false;
        }
    }

  return true;
}

SecondarySpecies::SecondarySpecies() {
    stoichiometry = {};
    equilibrium_constant = 1.0;
}

bool SecondarySpecies::operator==(const SecondarySpecies& other) {
    return (compare_doubles(this->equilibrium_constant, other.equilibrium_constant))
    && (this->stoichiometry == other.stoichiometry);
}

MineralSpecies::MineralSpecies() {
    stoichiometry = {};
    equilibrium_constant = 1.0;
    rate_constant = 1.0;
    molar_mass = 0.0;
    molar_volume = 0.0;
}

bool MineralSpecies::operator==(const MineralSpecies& other) {
    return (this->stoichiometry == other.stoichiometry) 
    && (compare_doubles(this->equilibrium_constant, other.equilibrium_constant))
    && (compare_doubles(this->molar_mass, other.molar_mass))
    && (compare_doubles(this->molar_volume, other.molar_volume))
    && (compare_doubles(this->rate_constant, other.rate_constant));
}

bool operator==(const map<string, SecondarySpecies>& l, const map<string, SecondarySpecies>& r) {
    if (l.size() != r.size()) return false;

    for (auto x : l) {
        if (!r.contains(x.first)) return false;

        if (r.at(x.first) != x.second) return false;
    }

    return true;    
}

bool operator==(const map<string, MineralSpecies>& l, const map<string, MineralSpecies>& r) {
    if (l.size() != r.size()) return false;

    for (auto x : l) {
        if (!r.contains(x.first)) return false;

        if (r.at(x.first) != x.second) return false;
    }

    return true;   
}

SecondarySpecies SecondarySpecies::from_yaml(const YAML::Node& node) {
    // Check for required sections
    if (!node[DBS_STOICHIOMETRY_TOKEN]) {
        std::cerr << "Stoichiometry section missing\n";
        exit(-1);
    }

    if (!node[DBS_LOGK_TOKEN]) {
        std::cerr << "Secondary species logK section missing\n";
        exit(-1);
    }
    
    const YAML::Node stoich_node = node[DBS_STOICHIOMETRY_TOKEN];
    double logK = node[DBS_LOGK_TOKEN].as<double>();

    map<string, double> stoich;
    for (YAML::const_iterator it = stoich_node.begin(); it != stoich_node.end(); it++) {
        stoich[it->first.as<string>()] = it->second.as<double>();
    }

    return SecondarySpecies(stoich, logK);
}

MineralSpecies MineralSpecies::from_yaml(const YAML::Node& node) {
    // Check for required species
    if (!node[DBS_STOICHIOMETRY_TOKEN]) {
        std::cerr << "Mineral species missing stoichiometry\n";
        exit(-1);
    }
    
    if (!node[DBS_LOGK_TOKEN]) {
        std::cerr << "Mineral species logK missing\n";
        exit(-1);
    }

    if (!node[DBS_MOLAR_MASS_TOKEN]) {
        std::cerr << "Mineral species molar mass entry missing\n";
        exit(-1);
    }

    if (!node[DBS_MOLAR_VOLUME_TOKEN]) {
        std::cerr << "Mineral species molar volume entry missing\n";
        exit(-1);
    }

    if (!node[DBS_RATE_CONST_TOKEN]) {
        std::cerr << "Mineral kinetic rate constant entry missing\n";
        exit(-1);
    }

    const YAML::Node stoich_node = node[DBS_STOICHIOMETRY_TOKEN];
    const double logK = node[DBS_LOGK_TOKEN].as<double>();
    const double rate_const = node[DBS_RATE_CONST_TOKEN].as<double>();
    const double molar_volume = node[DBS_MOLAR_VOLUME_TOKEN].as<double>();
    const double molar_mass = node[DBS_MOLAR_MASS_TOKEN].as<double>();

    map<string, double> stoich;
    for (YAML::const_iterator it = stoich_node.begin(); it != stoich_node.end(); it++) {
        stoich[it->first.as<string>()] = it->second.as<double>();
    }

    return MineralSpecies(stoich, molar_volume, molar_mass, rate_const, logK);
}

Database Database::from_string(const string& yamlString) {
    const YAML::Node node = YAML::Load(yamlString);

    // Check for sections
    if (!node[DBS_PRIMARY_SECTION_TOKEN]) {
        std::cerr << "Database missing primary species section\n";
        exit(-1);
    }

    if (!node[DBS_SECONDARY_SECTION_TOKEN]) {
        std::cerr << "Database missing secondary species section\n";
        exit(-1);
    }

    if (!node[DBS_MINERAL_SECTION_TOKEN]) {
        std::cerr << "Database missing mineral kinetic section\n";
        exit(-1);
    }

    const YAML::Node prim_node = node[DBS_PRIMARY_SECTION_TOKEN];
    const YAML::Node sec_node = node[DBS_SECONDARY_SECTION_TOKEN];
    const YAML::Node min_node = node[DBS_MINERAL_SECTION_TOKEN];

    // Produce each section
    vector<string> primary_species;
    map<string, SecondarySpecies> sec_species;
    map<string, MineralSpecies> min_species;

    for (YAML::const_iterator it = prim_node.begin(); it != prim_node.end(); it++) {
        primary_species.push_back(it->as<string>());
    }

    for (YAML::const_iterator it = sec_node.begin(); it != sec_node.end(); it++) {
        sec_species[it->first.as<string>()] = SecondarySpecies::from_yaml(it->second);
    }

    for (YAML::const_iterator it = min_node.begin(); it != min_node.end(); it++) {
        min_species[it->first.as<string>()] = MineralSpecies::from_yaml(it->second);
    }

    return Database(primary_species, sec_species, min_species);
}

vector<string> Database::get_primary_species() {
    return _primary_species;
}

map<string, SecondarySpecies> Database::get_secondary_species() {
    return _secondary_species;
}

map<string, MineralSpecies> Database::get_mineral_species() {
    return _mineral_kinetic;
}