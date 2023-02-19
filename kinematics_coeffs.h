#ifndef KINEMATICS_COEFFS_H
#define KINEMATICS_COEFFS_H

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

#include <boost/math/special_functions/gamma.hpp>

using namespace boost::math;

const double pi = 3.141592653589793238462643383279502884;


//////////////////////////////////////////////////////////////////////////////
// FROM ROY, SINHA, SURYA (2013). Discrete geometry of a small causal diamond.
// Helpers - Intermediate Coefficients

/**
 * @brief Gets the constant chi_k, as of Roy, Sinha, Surya 2013.
 * 
 * @param d - int, number of dimensions of the manifold 
 * @param k - int
 * @return double: Value of the constant
 */
inline
double chi_k(double d, double k)
{
    return 1/k * pow( tgamma(d+1)/2 , k-1) * tgamma(d/2)* tgamma(d)
    / ( tgamma(k*d/2) * tgamma((k+1)*d/2) );
}



/**
 * @brief Gets the constant \xi_0
 * 
 * @param d - double. Dimension of the manifold. 3 or 4.
 * @return double: coefficient value
 */
inline
double xi_0(double d)
{
    double vol;
    if ((int)d==3){
        vol = 2.*pi;
    }
    else if ((int)d==4){
        vol = 4.0*pi;
    }
    else{
        std::cout << "Dimension must be 3 or 4!\n";
    } 
    return vol / (d*(d-1)*std::pow(2,d-1)) ;
}



/**
 * @brief Calculates Q_k constant value
 * 
 * @param k Size of chains
 * @param d Dimensions of the manifold
 * @param C_k Number of chains of size k
 * @param rho Number density of causet. Default 1 means discrete units.
 * @return double
 */
inline
double Q_k(double k, double d, double C_k, double rho = 1.)
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
    return ((k+1)*d + 2)*Q_k;
}


/**
 * @brief Calculates value of K_k without knowledge of Qk
 * 
 * @param k - Chain size 
 * @param d - Dimensions of the manifold
 * @param C_k Number of chains of size k
 * @param rho Number density of causet. 1 means discrete units.
 * @return double
 */
inline
double K_k(double k, double d, double C_k, double rho)
{
    return ((k+1)*d + 2)*Q_k(k, d, C_k, rho);
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
 * @param rho - Density of causet sprinkling. 1 means discrete units.
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
#endif /*KINEMATICS_COEFFS_H */