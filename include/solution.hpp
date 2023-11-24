#pragma once

/** @file */

#include <tuple>
#include "armadillo"

using std::tuple;
using arma::Col;
using arma::Mat;
using arma::vec;
using arma::mat;

/**
 * @brief Data structure for storing the chemical state of the system
 * 
 */
struct ChemicalState {
    vec concentration;          /** The vector of concentrations of each species */
    vec total_concentration;    /** The vector of total concentrations of each species for mass conservation */

    /**
     * @brief Construct a new Chemical State object
     * 
     */
    ChemicalState() { }
    
    /**
     * @brief Construct a new Chemical State object
     * 
     * @param conc 
     * @param tot_conc 
     */
    ChemicalState(const vec& conc, const vec& tot_conc) : concentration(conc), total_concentration(tot_conc) { }
};

/**
 * @brief Data structure containing the constants relating to the equilibrium parameters
 */
struct EquilibriumConstants {
    mat stoichiometry_matrix;       /** Stoichiometry matrix with shape (number of species)x(number of secondary species)*/
    mat stoich_null_space;          /** Null space of the stoichiometry matrix, used in solving the system */
    vec conc_particular_solution;   /** Particular solution of the underdetermined equilibrium constants */
    vec log_eq_consts;              /** Vector of the Base-10 logarithms of the equilibrium constants */

    /**
     * @brief Construct a new Equilibrium Constants object
     */
    EquilibriumConstants() { }

    /**
     * @brief Construct a new Equilibrium Constants object
     * @param stoich_mat 
     * @param eq_consts 
     */
    EquilibriumConstants(const mat& stoich_mat, const vec& eq_consts);
};

/**
 * @brief Data structure containing the mass conservation matrix
 */
struct TotalConstants {
    mat tot_mat;     /** Matrix of shape (number of primary species)x(number of species) */

    /**
     * @brief Construct a new Total Constants object
     * 
     */
    TotalConstants() { }

    /**
     * @brief Construct a new Total Constants object
     * 
     * @param tot_mat 
     */
    TotalConstants(const mat& tot_mat) : tot_mat(tot_mat) { };
};

/**
 * @brief Data structure containing the kinetic reaction coefficients
 * 
 */
struct KineticConstants {
    mat kin_mat;       /** Kinetic stoichiometry matrix where each row is the formula for the ionic activity product for each mineral */
    vec kin_const;     /** Vector of the Base-10 logarithm of the kinetic rate constants */
    vec eq_const;      /** Vector of the base-10 logarithm of the equilibrium constants */

    /**
     * @brief Construct a new Kinetic Constants object
     * 
     */
    KineticConstants() { }

    /**
     * @brief Construct a new Kinetic Constants object
     * 
     * @param kin_mat 
     * @param kin_consts 
     * @param eq_const 
     */
    KineticConstants(const mat& kin_mat, const vec& kin_consts, const vec& eq_const) :
        kin_mat(kin_mat),
        kin_const(kin_consts),
        eq_const(eq_const) { }
};
