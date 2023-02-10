#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <set>
#include <stack>
#include <stdio.h>
#include <stdexcept>
#include <string>
#include <vector>
#include <chrono>
#include <unordered_set>
#include <chrono>

#include "../causets_cpp/sprinkledcauset.h"
#include "../causets_cpp/shapes.h"
#include "../causets_cpp/spacetimes.h"

#include "../causets_cpp/functions.h"
#include "../causets_cpp/vecfunctions.h"

#include <boost/range/combine.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <omp.h>

// $HOME var get
#include <unistd.h>
#include <sys/types.h>
//#include <pwd.h>

using namespace std::chrono;
using namespace boost::math;

const double pi = 3.1415926535897;



// Defining constants (functions) to be used later in the calculations for R...
///////////////////////////////////////////////////////////////////////////////

/**
 * @brief Gets the constant \xi_0
 * 
 * @param d - dimensions of the manifold 
 * @param r - radius of the (d-2)sphere
 * @return double: coefficient value
 */
double xi_0(double d, double r)
{
    
    double vol;
    if ((int)d==3){
        vol = 1;
    }
    else if (d==4){
        vol = pi*r*r;
    }
    else{
        std::cout << "Dimension must be 2,3 or 4!\n";
    } 
    return vol/(d*(d-1)*std::pow(2,d-1));
}


/**
 * @brief Gets the constant chi_k
 * 
 * @param d - number of dimensions of the manifold 
 * @param k - integer
 * @return double: Value of the constant
 */
double chi_k(double d, double k){
    return 1/k * pow(tgamma(d+1)/2, k-1) * tgamma(d/2)* tgamma(d)
    / ( tgamma(k*d/2) * tgamma((k+1)*d/2) );
}

/**
 * @brief Calculates Q_k constant value
 * 
 * @param k Size of chains
 * @param d Dimensions of the manifold
 * @param C_k Number of chains of size k between two points 
 * @param rho Number density of causet
 * @param r Radius of the d-2 dimensional sphere?
 * @return double
 */
double Q_k(double k, double d, double C_k, double rho, double r)
{
    return 1/std::pow(xi_0(d,r)*rho,3) * std::pow(C_k/chi_k(d,k), 3/k);
}


/**
 * @brief Calculates value of K_k
 * 
 * @param k - Chain size 
 * @param d - Dimensions of the manifold
 * @param Q_k - Constant value defined above
 * @return double 
 */
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
 * @param r Radius of the d-2 dimensional sphere?
 * @return double
 * @return double 
 */
double K_k(double k, double d, double C_k, double rho, double r)
{
    return ((k+1)*d+2)*Q_k(k,d,C_k,rho,r);
}


/**
 * @brief Calculates value of J_k
 * 
 * @param k - Chain size
 * @param d - Dimensions of the manifold
 * @param K_k - Constant value defined above
 * @return double 
 */
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
 * @param r - Radius of d-2 dim sphere?
 * @return double 
 */
double J_k(double k, double d, double C_k, double rho, double r)
{
    return K_k(k, d, Q_k(k, d, C_k, rho, r))*(k*d+2);
}


// Defining Ricci scalar (components) and proper time extension of the interval
///////////////////////////////////////////////////////////////////////////////


/**
 * @brief Proper time extension of the interval 
 * 
 * @param d - Manifold dim
 * @param C_k - Number of chains of size k
 * @param rho - Density of causet sprinkling
 * @param r - Radius of d-2 dim sphere?
 * @return double 
 */
double T(double d, double C_k, double rho, double r)
{
    double J1 = J_k(1,d, C_k, rho, r);
    double J2 = J_k(2,d, C_k, rho, r);
    double J3 = J_k(3,d, C_k, rho, r);

    return std::pow(1/(d*d)*(J1-2*J2+J3), 1/d);
}

/**
 * @brief Ricci scalar value at the centre of the interval
 * 
 * @param d - Manifold dim
 * @param C_k - Number of chains of size k
 * @param rho - Density of causet sprinkling
 * @param r - Radius of d-2 dim sphere?
 * @return double 
 */
double R(double d, double C_k, double rho, double r) 
{
    double K1 = K_k(1,d,C_k,rho,r);
    double K2 = K_k(2,d,C_k,rho,r);
    double K3 = K_k(3,d,C_k,rho,r);

    double J1 = J_k(1,d,K1);
    double J2 = J_k(2,d,K2);
    double J3 = J_k(3,d,K3);

    return -4*(d+2)*(2*d+2)*(3*d+2) *
        std::pow(2,2/(3*d)) *
        std::pow(d,(4/(3*d)-1)) * 
        (K1-2*K2+K3) / std::pow(J1-2*J2+J3, 1+2/(3*d));
}

/**
 * @brief (0,0) component of the Ricci tensorat the centre of the interval
 * 
 * @param d - Manifold dim
 * @param C_k - Number of chains of size k
 * @param rho - Density of causet sprinkling
 * @param r - Radius of d-2 dim sphere?
 * @return double 
 */
double R_00(double d, double C_k, double rho, double r) 
{
    double T_proper = T(d,C_k,rho,r);
    double Q1 = Q_k(1,d,C_k,rho,r);
    double Q2 = Q_k(2,d,C_k,rho,r);
    double Q3 = Q_k(3,d,C_k,rho,r);

    return -4*(2*d+2)*(3*d+2) / (std::pow(d,3) * std::pow(T_proper,3*d+2)) *
    ( (d+2)*Q1 - (5*d+4)*Q2 + (4*d+2)*Q3 );
}




int main(){

double x = xi_0(4,2);
double y = chi_k(4,2);
std::cout << "xi=" << x << std::endl;
std::cout << "chi=" << y << std::endl;

// end of main   
}