#include "potions.hpp"
#include "nlopt.h"

arma::mat EquilibriumSolution::stoichMat() {
    return arma::mat(this->_stoichMat);
}

arma::mat EquilibriumSolution::totConcMat() {
    return arma::mat(this->_totConcMat);
}

arma::vec EquilibriumSolution::eqConst() {
    return arma::vec(this->_equilibriumConstant);
}

arma::vec EquilibriumSolution::totConc() {
    return arma::vec(this->_totConc);
}

arma::vec EquilibriumSolution::primConc() {
    return arma::vec(this->_primConc);
}

arma::vec EquilibriumSolution::solve() {
    /* Solve for the equilibrium solution */
    arma::vec log_x = arma::log10(this->primConc());
    arma::mat null_stoich = arma::null(this->stoichMat());
    arma::vec log10_c_part = arma::pinv(this->stoichMat()) * log_x;    // Particular solution for the system
    auto log10_c = [&] (arma::vec log_x) {
        return log10_c_part + null_stoich * log_x;
    };

    auto get_c = [&] (arma::vec log_x) {
        return arma::pow(log10_c(log_x), 10.0);
    };

    // Define the residual function
    auto residual = [&] (arma::vec log_x) {
        return this->totConcMat() * get_c(log_x) - this->totConc();
    };

    // Find the root of the matrix
    throw NotImplemented();
}