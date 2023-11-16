#include <string>
#include <cmath>
#include "catch2/catch_test_macros.hpp"
#include "armadillo"
#include "potions.hpp"

using std::string;

// Static strings
const string CHEMFILE_STR = R"(SimulationType:
  Equilibrium   # Other option is "Kinetic"
PrimarySpecies:
  H+: -999      # The value H+ is arbitrary
  Ca++: 1e-3
  HCO3-: 2e-3
SecondarySpecies:
  - CO3--
  - CO2(aq)
MineralSpecies:
  Calcite: 1.0  # Surface area: m^2/g)";

const string DATABASE_STR = R"(
PrimarySpecies:
  - Ca++
  - H+
  - HCO3-
SecondarySpecies:
  CO3--:
    stoichiometry:
      H+: 1
      HCO3-: -1
    logK: 9.617
  CO2(aq):
    stoichiometry:
      H+: -1
      HCO3-: -1
    logK: -6.3447
MineralSpecies:
  Calcite:
    stoichiometry:
      Ca++: 1
      HCO3-: 1
      H+: -1
    logK: -7.30
    rate: -9.19
    molarVolume: 2.5
    molarMass: 2.5)";

const string secString1 = R"(stoichiometry:
  H+: -1
  HCO3-: 1
logK: 9.617)";

const string secString2 = R"(CO2(aq):
stoichiometry:
    H+: 1
    HCO3-: 1
logK: -6.3447)";

const string minString = R"(CaCO3:
stoichiometry:
  Ca++: 1
  HCO3-: 1
  H+: -1
logK: -7.30
rate: -9.19
molarVolume: 2.5
molarMass: 2.5)";

bool compare_float(double a, double b) {
    return std::abs(a - b) < 1e-12;
}

bool compare_maps(const map<string,double>& a, const map<string,double>& b) {
    if (a.size() != b.size()) return false;

    for (auto a_i: a) {
        if (!b.contains(a_i.first)) {
            return false;
        }

        if (!compare_float(a_i.second, b.at(a_i.first))) {
        return false;
        }
    }

  return true;
}

TEST_CASE("Test reading 'chem.yaml'", "[ChemFile::FromFile]") {
    ChemFile chem = ChemFile::FromString(CHEMFILE_STR);

    const map<string, double> CHEM_PRIM_SPECIES = {
      {"H+", -999.0 }, 
      {"Ca++", 1e-3}, 
      {"HCO3-", 2e-3}
      };
    const vector<string> CHEM_SEC_SPECIES = {"CO3--", "CO2(aq)"};
    const map<string, double> CHEM_MIN_SPECIES = {{"Calcite", 1.0}};

    REQUIRE(chem.runType() == PotionsRunType::EQUILIBRIUM);

    REQUIRE(compare_maps(chem.primarySpecies(), CHEM_PRIM_SPECIES));
    REQUIRE(chem.secondarySpecies() == CHEM_SEC_SPECIES);
    REQUIRE(compare_maps(chem.mineralSpecies(), CHEM_MIN_SPECIES));
}

TEST_CASE("[SecondarySpecies::FromYaml]") {
    {
        // Case 1
        const YAML::Node node1 = YAML::Load(secString1);
        map<string, double> s1_stoich = {{"H+",-1.0}, {"HCO3-", 1.0}};
        SecondarySpecies s1 = SecondarySpecies::FromYaml(node1);
        REQUIRE(s1.stoichiometry == s1_stoich);
        REQUIRE(compare_float(s1.equilibriumConstant, 9.617));
    }

    {
        // Case 2
        YAML::Node node2 = YAML::Load(secString2);
        map<string, double> s2_stoich = {{"H+", 1.0}, {"HCO3-", 1.0}};
        SecondarySpecies s2 = SecondarySpecies::FromYaml(node2);
        REQUIRE(s2.stoichiometry == s2_stoich);
        REQUIRE(compare_float(s2.equilibriumConstant, -6.3447));
    }

}

TEST_CASE("[MineralSpecies::FromYaml]") {
    {
        const YAML::Node node = YAML::Load(minString);
        map<string,double> stoich = {
            {"Ca++", 1},
            {"HCO3-", 1},
            {"H+", -1}
        };

        MineralSpecies minSpec = MineralSpecies::FromYaml(node);

        REQUIRE(compare_maps(minSpec.stoichiometry, stoich));
        REQUIRE(compare_float(minSpec.rateConstant, -9.19));
        REQUIRE(compare_float(minSpec.equilibriumConstant, -7.30));
    }
}

TEST_CASE("Test reading 'cdbs.yaml'", "[Database::FromString]") {
    
    Database cdbs = Database::FromString(DATABASE_STR);

    // Define expected values
    vector<string> expectedPrimSpec = {"Ca++", "H+", "HCO3-"};
    map<string, double> stoich_1 = {{"H+", 1}, {"HCO3-", -1.0}};
    map<string, double> stoich_2 = {{"H+", -1.0}, {"HCO3-", -1.0}};
    map<string, SecondarySpecies> expectedSecSpec = {
      {"CO3--", SecondarySpecies(stoich_1, 9.617)},
      {"CO2(aq)", SecondarySpecies(stoich_2, -6.3447)}
      };
    map<string, MineralSpecies> expectedMinSpec = {
      {"Calcite", MineralSpecies({{"Ca++", 1}, {"HCO3-", 1}, {"H+", -1.0}}, 2.5, 2.5, -9.19, -7.30)}
    };

    REQUIRE(cdbs.getPrimarySpecies() == expectedPrimSpec);
    REQUIRE(cdbs.getSecondarySpecies() == expectedSecSpec);
    REQUIRE(cdbs.getMineralSpecies() == expectedMinSpec);
}


// Model inputs
TEST_CASE("[ModelInputs::initChemState]") {
    ChemFile chem = ChemFile::FromString(CHEMFILE_STR);
    Database cdbs = Database::FromString(DATABASE_STR);
    ModelInputs mod = ModelInputs(chem, cdbs);

    vec tot_conc = {1e-3, 0.0, 2e-3};

    ChemicalState chms = mod.initChemState();
    
    REQUIRE(chms.totalConcentration.size() == 3);

    vec err = tot_conc - chms.totalConcentration;
    REQUIRE(arma::abs(err).max() < 1e-12);
    REQUIRE(arma::abs(chms.concentration).max() < 1e-18);
    REQUIRE(chms.concentration.size() == 5);
}

TEST_CASE("[ModelInputs::equilibriumConstants]") {
    // Species map: H+, Ca++, HCO3-, CO3--, CO2(aq)
    const arma::mat stoich_mat = {
        {0, 1, -1, 1, 0},   // CO3--
        {0, -1, -1, 0, 1}   // CO2(aq)
    };
    const arma::vec log_k = {9.617, -6.3447};

    ChemFile chem = ChemFile::FromString(CHEMFILE_STR);
    Database cdbs = Database::FromString(DATABASE_STR);
    ModelInputs mod = ModelInputs(chem, cdbs);

    EquilibriumConstants eq = mod.equilibriumConstants();

    REQUIRE(eq.stoich_mat.n_cols == 5);
    REQUIRE(eq.stoich_mat.n_rows == 2);
    REQUIRE(eq.eq_consts.size() == 2);

    const arma::mat err_stoich = stoich_mat - eq.stoich_mat;
    const arma::vec err_log_k = log_k - eq.eq_consts;

    REQUIRE(arma::abs(err_stoich).max() < 1e-12);
    REQUIRE(arma::abs(err_log_k).max() < 1e-12);
}

TEST_CASE("[ModelInputs::totalConstants]") {
    // Species map: Ca++, H+, HCO3-, CO3--, CO2(aq)
    const arma::mat tot_mat = {
        {1, 0, 0, 0, 0},    // Ca++
        {2, 1, -1, -2, 0},  // H+ (charge balance)
        {0, 0, 1, 1, 1}     // HCO3-
    };

    ChemFile chem = ChemFile::FromString(CHEMFILE_STR);
    Database cdbs = Database::FromString(DATABASE_STR);
    ModelInputs mod = ModelInputs(chem, cdbs);

    const arma::mat calc_mat = mod.totalConstants().tot_mat;
    const arma::mat err = tot_mat - calc_mat;

    REQUIRE(calc_mat.n_rows == 3);
    REQUIRE(calc_mat.n_cols == 5);
    REQUIRE(arma::abs(err).max() < 1e-12);
}

TEST_CASE("[ModelInputs::kineticConstants]") {
    // Species map: Ca++, H+, HCO3-, CO3--, CO2(aq)
    const arma::mat kin_mat = {
        {1, -1, 1, 0, 0}
    };
    const arma::vec eq_consts = {-7.30};
    const arma::vec kin_consts = {-9.19};
    
    ChemFile chem = ChemFile::FromString(CHEMFILE_STR);
    Database cdbs = Database::FromString(DATABASE_STR);
    ModelInputs mod = ModelInputs(chem, cdbs);

    KineticConstants kin = mod.kineticConstants();

    REQUIRE(kin.kin_mat.n_rows == 1);
    REQUIRE(kin.kin_mat.n_cols == 5);
    REQUIRE(kin.eq_const.size() == 1);
    REQUIRE(kin.kin_const.size() == 1);

    REQUIRE(compare_float(eq_consts[0], kin.eq_const[0]));
    REQUIRE(compare_float(kin_consts[0], kin.kin_const[0]));

    const arma::mat err = kin_mat - kin.kin_mat;
    REQUIRE(arma::abs(err).max() < 1e-12);
}