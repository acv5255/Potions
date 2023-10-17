#pragma once

#include "armadillo"

class Solution {
    public:
        virtual arma::vec solve() = 0;
        virtual arma::mat stoichMat() = 0;
        virtual arma::vec primConc() = 0;
};

class EquilibriumSolution : public Solution {
    private:
        arma::mat _stoichMat;           // Stoichiometry matrix
        arma::mat _totConcMat;          // Total concentration matrix
        arma::vec _equilibriumConstant; // Equilibrium constant vector
        arma::vec _totConc;             // Total concentration in our system
        arma::vec _primConc;         // Primary concentration vector

    public:
        arma::vec solve() override;
        arma::mat stoichMat() override;
        arma::mat totConcMat();
        arma::vec eqConst();
        arma::vec totConc();
        arma::vec primConc() override;
};