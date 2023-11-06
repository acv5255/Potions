#include "kinetic.hpp"


template<typename T>
Col<T> kinetic_rate(const Col<T>& conc, const Col<T> surface_area, const KineticConstants<T>& kin_params) {
    const Col<T> log_c = arma::log(conc);
    const Col<T> log_iap = kin_params.kin_mat * log_c;
    const Col<T> iap = arma::exp(log_iap);
    return surface_area * kin_params.rate_consts * (1.0 - iap / kin_params.eq_consts);
}