#include "kinetic.hpp"
#include "ode.hpp"


vec kinetic_rate(const vec& totConc, const vec& surface_area, const EquilibriumConstants& eq, const KineticConstants& kin, const TotalConstants& tot) {
    /*
        Takes in concentrations and returns the rate of the change of hte total concentrations
     */
    ChemicalState chem = solve_equilibrium(totConc, eq, tot);

    const vec log_c = arma::log(chem.concentration);

    const vec log_iap = kin.kin_mat * log_c;
    vec kin_consts_corrected = arma::zeros(kin.kin_const.size());
    for (int i = 0; i < kin.kin_const.size(); i++) {
        kin_consts_corrected[i] = std::pow(10.0, kin.kin_const[i]);
    }
    vec eq_consts_corrected = arma::zeros(kin.eq_const.size());
    for (int i = 0; i < kin.eq_const.size(); i++) {
        eq_consts_corrected[i] = std::pow(10.0, kin.eq_const[i]);
    }

    const vec iap = arma::exp(log_iap);
    const vec mineralRates = surface_area * kin_consts_corrected * (1.0 - iap / eq_consts_corrected);
    const mat minRateMat = arma::diagmat(mineralRates);
    const mat componentMat = minRateMat * kin.kin_mat;
    const vec speciesRate = componentMat.t() * arma::ones(kin.kin_mat.n_rows);

    return tot.tot_mat * speciesRate;
}

ChemicalState solve_kinetic_equilibrium(
    const ChemicalState& chem,
    const vec& surfaceArea,
    const KineticConstants& kin,
    const EquilibriumConstants& eq,
    const TotalConstants& tot,
    double dt
) {
    // Construct the ODE
    auto getKineticRate = [&] (const vec& totConc) {
        return kinetic_rate(totConc, surfaceArea, eq, kin, tot);
    };

    auto func_val = getKineticRate(chem.total_concentration);

    vec new_tot_conc = solve_ode<double>(getKineticRate, chem.total_concentration, dt);

    return solve_equilibrium(new_tot_conc, eq, tot);
}