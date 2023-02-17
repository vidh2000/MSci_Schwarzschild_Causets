#ifndef KINEMATICS_FUNCTIONS_H
#define KINEMATICS_FUNCTIONS_H

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

#include "causet.h"
#include "embeddedcauset.h"
#include "kinematics_coeffs.h"
#include "functions.h"
#include "vecfunctions.h"

#include <boost/range/combine.hpp>
#include <boost/math/special_functions/gamma.hpp>

using namespace boost::math;


/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
// Causet Sampling Functions for Curvature Information
/////////////////////////////////////////////////////////////////////////////

/**
 * @brief Sample Ricci scalar and R00 values in DISCRETEENSS UNITS 
 * Nsamples times in the causet 
 * by using the estimator from RSS - Roy, Sinha, Surya 2013 paper - 
 * at the centre of various intervals. 
 * It assumes that coord[1] is some radial distance.
 * 
 * @param Causet EmbeddedCauset : Causet to be sampled.
 * @param Nsamples int : number of samples to analyse.
 * @param interval_sizemin int : minimum size of intervals to take into
 * consideration.
 * @param interval_sizemax int : maximum size of intervals to take into 
 * account. Default is 0, meaning the causet's size.
 * 
 * @exception std::runtime_error - if does not find suitable interval in N_max 
 * tries.
 * 
 * @return std::vector<std::vector <double> > (Nsamples, 3): vector of 3-vectors
 * with Ricci Scalar, Ricci 00 (as of RSS 2013) and radial distance.
 */
inline
std::vector<std::vector <double>> Riccis_RSS_radial_sample
(EmbeddedCauset & Causet, int Nsamples, 
int interval_sizemin, int interval_sizemax = 0) 
{
    std::vector<std::vector <double>> vector_of_Ricci_Ricci00_radius
    (Nsamples);

    int d = Causet._spacetime._dim;
    double rho = 1.;

    std::vector< std::pair< std::vector<double>, double> > 
    intervals_chains_and_r = 
    Causet.get_Nchains_inInterval(Nsamples, interval_sizemin, 
                                    4, interval_sizemax);
    
    #pragma omp paraller for
    for (int n = 0; n<Nsamples; n++)
    {
        std::vector<double> C_k = intervals_chains_and_r[n].first;

        double Q1 = Q_k(1,d,C_k[1],rho);
        double Q2 = Q_k(2,d,C_k[2],rho);
        double Q3 = Q_k(3,d,C_k[3],rho);

        double K1 = K_k(1,d,Q1);
        double K2 = K_k(2,d,Q2);
        double K3 = K_k(3,d,Q3);

        double J1 = J_k(1,d,K1);
        double J2 = J_k(2,d,K2);
        double J3 = J_k(3,d,K3);

        double T_proper = std::pow(1/(d*d)*(J1-2*J2+J3), 1/(3*d));
        
        double R00_n =-4*(2*d+2)*(3*d+2) /
        (std::pow(d,3) * std::pow(T_proper,3*d+2)) *
        ( (d+2)*Q1 - (5*d+4)*Q2 + (4*d+2)*Q3 );

        double R_RSS_n = -4*(d+2)*(2*d+2)*(3*d+2) *
                            std::pow(2,2/(3*d)) *
                            std::pow(d,(4/(3*d)-1)) * 
                            (K1-2*K2+K3) / std::pow(J1-2*J2+J3, 1+2/(3*d));
        
        vector_of_Ricci_Ricci00_radius[n] =
                            {R_RSS_n, R00_n, intervals_chains_and_r[n].second};   
    }
    return vector_of_Ricci_Ricci00_radius;
}


/**
 * @brief Sample Ricci scalar value Nsamples times in the causet by using the 
 * estimator from RSS - Roy, Sinha, Surya 2013 paper - 
 * at the centre of various intervals. 
 * It assumes that coord[1] is some radial distance.
 * 
 * @param Causet EmbeddedCauset : Causet to be sampled.
 * @param Nsamples int : number of samples to analyse.
 * @param interval_sizemin int : minimum size of intervals to take into
 * consideration.
 * @param interval_sizemax int : maximum size of intervals to take into 
 * account. Default is 0, meaning the causet's size.
 * 
 * @exception std::runtime_error - if does not find suitable interval in N_max 
 * tries.
 * 
 * @return std::vector<std::pair <double, double> > : vector of pairs 
 * Ricci Scalar (as of RSS 2013) and radial distance 
 */
inline
std::vector<std::pair <double, double> > R_RSS_radial_sample
(EmbeddedCauset & Causet, int Nsamples, 
int interval_sizemin, int interval_sizemax = 0) 
{
    std::vector <std::pair <double, double > >
    vector_of_pairs_of_RicciScalar_and_radius = {};

    int d = Causet._spacetime._dim;

    std::vector< std::pair< std::vector<double>, double> > 
    intervals_chains_and_r = 
    Causet.get_Nchains_inInterval(Nsamples, interval_sizemin, 
                                    4, interval_sizemax);
    for (int n = 0; n<Nsamples; n++)
    {
        std::vector<double> C_k = intervals_chains_and_r[n].first;

        double K1 = K_k(1,d,C_k[0],1.);
        double K2 = K_k(2,d,C_k[1],1.);
        double K3 = K_k(3,d,C_k[2],1.);

        double J1 = J_k(1,d,K1);
        double J2 = J_k(2,d,K2);
        double J3 = J_k(3,d,K3);

        double R_RSS_n = -4*(d+2)*(2*d+2)*(3*d+2) *
                            std::pow(2,2/(3*d)) *
                            std::pow(d,(4/(3*d)-1)) * 
                            (K1-2*K2+K3) / std::pow(J1-2*J2+J3, 1+2/(3*d));
        
        vector_of_pairs_of_RicciScalar_and_radius.push_back(
                                {R_RSS_n, intervals_chains_and_r[n].second});   
    }
    
    return vector_of_pairs_of_RicciScalar_and_radius;
}


/**
 * @brief Sample Ricci Tensor 00 value Nsample times in the causet by using the 
 * estimator from RSS - Roy, Sinha, Surya 2013 paper - 
 * at the centre of various intervals. 
 * It assumes that coord[1] is some radial distance.
 * 
 * @param Causet EmbeddedCauset : Causet to be sampled.
 * @param Nsamples int : number of samples to analyse.
 * @param interval_sizemin int : minimum size of intervals to take into
 * consideration.
 * @param interval_sizemax int : maximum size of intervals to take into 
 * account. Default is 0, meaning the causet's size.
 * 
 * @exception std::runtime_error - if does not find suitable interval in N_max 
 * tries.
 * 
 * @return std::vector<std::pair <double, double> > : vector of pairs 
 * Ricci Scalar (as of RSS 2013) and radial distance 
 */
inline
std::vector<std::pair <double, double> > R00_RSS_radial_sample
(EmbeddedCauset & Causet, int Nsamples, 
int interval_sizemin, int interval_sizemax = 0) 
{
    std::vector <std::pair <double, double > >
    vector_of_pairs_of_Ricci00_and_radius = {};

    int d = Causet._spacetime._dim;
    double rho = 1.;  // in discreteness units

    std::vector< std::pair< std::vector<double>, double> > 
    intervals_chains_and_r = 
    Causet.get_Nchains_inInterval(Nsamples, interval_sizemin, 
                                    4, interval_sizemax);
    for (int n = 0; n<Nsamples; n++)
    {
        std::vector<double> C_k = intervals_chains_and_r[n].first;
        
        double Q1 = Q_k(1,d,C_k[1],rho);
        double Q2 = Q_k(2,d,C_k[2],rho);
        double Q3 = Q_k(3,d,C_k[3],rho);

        double K1 = K_k(1,d,Q1);
        double K2 = K_k(2,d,Q2);
        double K3 = K_k(3,d,Q3);

        double J1 = J_k(1,d,K1);
        double J2 = J_k(2,d,K2);
        double J3 = J_k(3,d,K3);

        double T_proper = std::pow(1/(d*d)*(J1-2*J2+J3), 1/(3*d));

        double R00_n =-4*(2*d+2)*(3*d+2) /
        (std::pow(d,3) * std::pow(T_proper,3*d+2)) *
        ( (d+2)*Q1 - (5*d+4)*Q2 + (4*d+2)*Q3 );
        
        vector_of_pairs_of_Ricci00_and_radius.push_back(
                                {R00_n, intervals_chains_and_r[n].second});   
    }
    
    return vector_of_pairs_of_Ricci00_and_radius;
}


/**
 * @brief Sample causet for Ricci scalar as from Benincasa-Dowker operator, 
 * in DISCRETENESS UNITS. 
 * Assumes natural labelling, with coords[1] being some sort of radial 
 * coordinate from origin and only samples in central two quartiles 
 * of time and space to avoid boundaries effects. 
 * 
 * @param Causet EmbeddedCauset to be sampled.
 * @param Nsamples int : number of samples.
 * @return std::vector<std::pair <double, double> > : vector of pairs 
 * {RicciScalar, radial distance}.
 */
inline
std::vector<std::pair <double, double> > R_BD_sample
(EmbeddedCauset & Causet, int Nsamples)
{
    std::vector<std::pair <double, double> > vec_of_pair_R_r(Nsamples);

    std::vector<double> center = Causet._shape._center; 
    double duration = Causet._shape._params.find("duration")->second;
    double radius   = Causet._shape._params.find( "radius" )->second;
    double tmin = center[0] - 0.25*duration;
    double tmax = center[0] + 0.25*duration;
    double rmin = 0.25*radius;
    double rmax = 0.75*radius;

    // Define mersenne_twister_engine Random Gen. (with random seed)
    std::random_device rd;
    int seed = rd();
    std::mt19937 gen(seed);
    std::uniform_real_distribution<> dis(0, Causet._size);

    #pragma omp paraller for
    for (int rep = 0; rep < Nsamples; rep ++)
    {
        int xi = (int) dis(gen);
        double ri = Causet._coords[xi][1];
        double ti = Causet._coords[xi][0];

        while (!(tmin <= ti && ti <= tmax && rmin <= ri && ri <= rmax)){
            xi += (int) dis(gen) - xi;
            ri += Causet._coords[xi][1] + ri ;
            ti += Causet._coords[xi][0] - ti;
        }
        
        std::vector<double> N_arr = Causet.Nk_BD(xi, 4, 1);
        double Ri = 4*std::pow(2./3.,0.5) *
                (1 - (N_arr[0] - 9*N_arr[1] + 16*N_arr[2] - 8*N_arr[3]) );
        
        vec_of_pair_R_r[rep] = {Ri*1., ri*1.};
    }

    return vec_of_pair_R_r;
}


//////////////////////////////////////////////////////////////////////////////
// FROM ROY, SINHA, SURYA (2013). Discrete geometry of a small causal diamond.
// Single Interval Functions
//////////////////////////////////////////////////////////////////////////////

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

    return std::pow(1/(d*d)*(J1-2*J2+J3), 1/(3*d));
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
    print(C_k_arr);
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
// Benincasa Dowker
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


/**
 * @brief Ricci scalar as calculated from Benincasa-Dowker operator. It is 
 * already in discreteness units. 
 * 
 * @param x int : Label of element at which compute BD Ricci Scalar. 
 * @param C Causet : Causet to which x belongs.
 * 
 * @return double 
 */
inline
double R_BD(int x, Causet & C)
{
    std::vector<double> N_arr = C.Nk_BD(x, 4, 1);
    return 4*std::pow(2./3.,0.5) *
            (1 - (N_arr[0] - 9*N_arr[1] + 16*N_arr[2] - 8*N_arr[3]) );
}


/**
 * @brief Ricci scalar as calculated from Benincasa-Dowker operator. It is
 * already in discreteness units.
 * 
 * @param x int : Label of element at which compute BD Ricci Scalar. 
 * @param N_arr vector<double> : {N1, N2, N3, N4} where
 *          N_k: Number of y in causet such that |I(y,x)|= k + 1
 * @return double 
 */
inline
double R_BD(int x,  std::vector<double> N_arr)
{
    return 4*std::pow(2,0.5)/(std::pow(3,0.5)) *
            (1-(N_arr[0]-9*N_arr[1]+16*N_arr[2]-8*N_arr[3]));
}



#endif /* KINEMATICS_FUNCTIONS_H */