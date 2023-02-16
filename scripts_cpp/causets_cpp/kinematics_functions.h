#ifndef SCHWARZSCHILD_KINEMATICS_H
#define SCHWARZSCHILD_KINEMATICS_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <fstream>
#include <stack>
#include <string>
#include <stdio.h>
#include <vector>
#include <set>
#include <unordered_set>

#include "functions.h"
#include "vecfunctions.h"

#include <boost/range/combine.hpp>
#include <boost/math/special_functions/gamma.hpp>

using namespace boost::math;

const double pi = 3.141592653589793238462643383279502884;
                  
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// Defining constants (functions) to be used later in various kin. calculations
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


/**
 * @brief Gets the constant chi_k
 * 
 * @param d - number of dimensions of the manifold 
 * @param k - integer
 * @return double: Value of the constant
 */
inline
double chi_k(double d, double k){
    return 1/k * pow(tgamma(d+1)/2, k-1) * tgamma(d/2)* tgamma(d)
    / ( tgamma(k*d/2) * tgamma((k+1)*d/2) );
}


/**
 * @brief Gets the constant \xi_0
 * 
 * @param d - dimensions of the manifold 
 * @return double: coefficient value
 */
inline
double xi_0(double d)
{
    
    double vol;
    if ((int)d==3){
        vol = pi;
    }
    else if (d==4){
        vol = 4.0/3.0*pi;
    }
    else{
        std::cout << "Dimension must be 3 or 4!\n";
    } 
    return vol/(d*(d-1)*std::pow(2,d-1));
}



/**
 * @brief Calculates Q_k constant value
 * 
 * @param k Size of chains
 * @param d Dimensions of the manifold
 * @param C_k Number of chains of size k between two points 
 * @param rho Number density of causet
 * @return double
 */
inline
double Q_k(double k, double d, double C_k, double rho)
{
    return 1/std::pow(xi_0(d)*rho,3) * std::pow(C_k/chi_k(d,k), 3/k);
}

/**
 * @brief Calculates value of K_k
 * 
 * @param k - Chain size 
 * @param d - Dimensions of the manifold
 * @param Q_k - Constant value defined above
 * @return double 
 */
inline
double K_k(double k, double d, double Q_k)
{
    return ((k+1)*d+2)*Q_k;
}


/**
 * @brief Overridden: Calculates value of K_k without knowledge of Q_k
 * 
 * @param k Size of chains
 * @param d Dimensions of the manifold
 * @param C_k Number of chains of size k between two points 
 * @param rho Number density of causet
 * @return double
 * @return double 
 */
inline
double K_k(double k, double d, double C_k, double rho)
{
    return ((k+1)*d+2)*Q_k(k,d,C_k,rho);
}


/**
 * @brief Calculates value of J_k
 * 
 * @param k - Chain size
 * @param d - Dimensions of the manifold
 * @param K_k - Constant value defined above
 * @return double 
 */
inline
double J_k(double k, double d, double K_k)
{
    return K_k*(k*d+2);
}


/**
 * @brief Overridden J_k calculating the value from without the knowledge
 *          of the values of other constants
 * 
 * @param k - Chain size
 * @param d - Manifold dim
 * @param C_k - Number of chains of size k
 * @param rho - Density of causet sprinkling
 * @return double 
 */
inline
double J_k(double k, double d, double C_k, double rho)
{
    return K_k(k, d, Q_k(k, d, C_k, rho))*(k*d+2);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// MM_dim estimator for curved spacetimes
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/**
 * @brief Calculate n choose k and return the value as double
 * 
 * @param n Integer
 * @param k Integer
 * @return double 
 */
inline
double binomialCoefficient(int n, int k)
{
    int C[n + 1][k + 1];
    int i, j;

    // Caculate value of Binomial Coefficient in bottom up manner
    for (i = 0; i <= n; i++) {
    for (j = 0; j <= std::min(i, k); j++) {
        // Base Cases
        if (j == 0 || j == i) {
        C[i][j] = 1;
        } else {
        // Calculate value using previously stored values
        C[i][j] = C[i - 1][j - 1] + C[i - 1][j];
        }
    }
    }
    return (double)C[n][k];
}



/**
 * @brief 
 * 
 * @param d - Manifold dimension 
 * @param C_k_arr - List of C_ks (Numbers of k-long chains)
 *               for chain lengths k=1,2,3,4.
 * @return double
 */
inline
double MMdim_eqn(double d, std::vector<double> C_k_arr)
{
    // C_k_arr must have length 4
    if (C_k_arr.size() !=4)
    {
        std::cout << "C_k_arr must contain #chains for k=1,2,3,4." << std::endl;
    }
    std::vector<double> k_arr = {1.0,2.0,3.0,4.0};
    double result = 0;
    
    for (auto && tup : boost::combine(k_arr, C_k_arr))
    {
        double k, C_k;
        boost::tie(k,C_k) = tup;

        result += std::pow(-1,k)*binomialCoefficient(3,k-1)*
            (k*d+2)*((k+1)*d+2)*std::pow(C_k, 4/k) /
            std::pow(chi_k(d,k), 4/k);
    }
    return result;
}


/**
 * @brief Yields the estimate of the MM-dimension for curved spacetime
 * 
 * @param d - Dimension of manifold
 * @param C_k_arr - array of chain lengths for 1,2,3,4-chains
 * @return double 
 */
inline
double estimate_MMd(std::vector<double> C_k_arr)
{
    // Define function whose root needs to be found
    auto MM_to_solve = [C_k_arr](double d){
        return MMdim_eqn(d,C_k_arr);
    };

    double dmin = 0.1;
    double dmax = 10;
    // Estimate dimension of Causet
    double dim_estimate = bisection(MM_to_solve,dmin,dmax);
    return dim_estimate;
};




///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// Functions for finding the Ricci scalar and tensor and proper time T
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////





/**
 * @brief Proper time extension of the interval 
 * 
 * @param d - Manifold dim
 * @param C_k - Vector of: Number of chains of size k for k=1,2,3
 * @param rho - Density of causet sprinkling
 * @return double 
 */
inline
double T(double d, std::vector<double> C_k, double rho)
{
    double J1 = J_k(1,d, C_k[0], rho);
    double J2 = J_k(2,d, C_k[1], rho);
    double J3 = J_k(3,d, C_k[2], rho);

    return std::pow(1/(d*d)*(J1-2*J2+J3), 1/d);
}


/**
 * @brief Ricci scalar value at the centre of the interval
 *          (from RSS - Roy, Sinha, Surya 2013 paper)
 * @param d - Manifold dim
 * @param C_k - Vector of: Number of chains of size k for k=1,2,3
 * @param rho - Density of causet sprinkling
 * @return double 
 */
inline
double R_RSS(double d, std::vector<double> C_k, double rho) 
{
    double K1 = K_k(1,d,C_k[0],rho);
    double K2 = K_k(2,d,C_k[1],rho);
    double K3 = K_k(3,d,C_k[2],rho);

    double J1 = J_k(1,d,K1);
    double J2 = J_k(2,d,K2);
    double J3 = J_k(3,d,K3);

    return -4*(d+2)*(2*d+2)*(3*d+2) *
        std::pow(2,2/(3*d)) *
        std::pow(d,(4/(3*d)-1)) * 
        (K1-2*K2+K3) / std::pow(J1-2*J2+J3, 1+2/(3*d));
}




/**
 * @brief (0,0) component of the Ricci tensor at the centre of the interval
 *              (from RSS - Roy, Sinha, Surya 2013 paper)
 * @param d - Manifold dim
 * @param C_k - Vector of: Number of chains of size k for k=1,2,3
 * @param rho - Density of causet sprinkling
 * @return double 
 */
inline
double R_00(double d, std::vector<double> C_k, double rho) 
{
    double T_proper = T(d,C_k,rho);
    print(T_proper);
    double Q1 = Q_k(1,d,C_k[1],rho);
    double Q2 = Q_k(2,d,C_k[2],rho);
    double Q3 = Q_k(3,d,C_k[3],rho);

    return -4*(2*d+2)*(3*d+2) / (std::pow(d,3) * std::pow(T_proper,3*d+2)) *
    ( (d+2)*Q1 - (5*d+4)*Q2 + (4*d+2)*Q3 );
}



/**
 * @brief Ricci scalar as calculated via Benincasa-Dowker action
 * 
 * @param l - discreteness length
 * @param N_arr - array (N1, N2, N3, N4) where
 *          N_i: Number of i-chains in the interval
 * @return double 
 */
inline
double R_BD(double l, std::vector<double> N_arr)
{
    return 4*std::pow(2,0.5)/(std::pow(3,0.5)*l*l) *
            (1-(N_arr[0]-9*N_arr[1]+16*N_arr[2]-8*N_arr[3]));
}



#endif /* SCHWARZSCHILD_KINEMATICS_H */