#include "kinetic.hpp"


vec kinetic_rate(const vec& conc, const vec& surface_area, const KineticConstants& kin_params) {
    const vec log_c = arma::log(conc);
    const vec log_iap = kin_params.kin_mat * log_c;
    const vec iap = arma::exp(log_iap);
    return surface_area * kin_params.kin_const * (1.0 - iap / kin_params.eq_const);
}