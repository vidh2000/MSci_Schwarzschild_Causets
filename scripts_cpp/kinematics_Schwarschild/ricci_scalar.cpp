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
 *          (from RSS - Roy, Sinha, Surya 2013 paper)
 * @param d - Manifold dim
 * @param C_k - Number of chains of size k
 * @param rho - Density of causet sprinkling
 * @param r - Radius of d-2 dim sphere?
 * @return double 
 */
double R_RSS(double d, double C_k, double rho, double r) 
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
 * @brief (0,0) component of the Ricci tensor at the centre of the interval
 *              (from RSS - Roy, Sinha, Surya 2013 paper)
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


/**
 * @brief Ricci scalar as calculated via Benincasa-Dowker action
 * 
 * @param l - discreteness length
 * @param N_arr - array (N1, N2, N3, N4) where
 *          N_i: Number of i-chains in the interval
 * @return double 
 */
double R_BD(double l, std::vector<double> N_arr)
{
    return 4*std::pow(2,0.5)/(std::pow(3,0.5)*l*l) *
            (1-(N_arr[0]-9*N_arr[1]+16*N_arr[2]-8*N_arr[3]));
}


int main(){

double x = xi_0(4,2);
double y = chi_k(4,2);
std::cout << "xi=" << x << std::endl;
std::cout << "chi=" << y << std::endl;


// Create causet in Schwarzschild spacetime
///////////////////////////////////////////////////////////////////////////////
                
std::vector<int> dims = {4}; 
std::vector<int> cards = {1000};
double mass = 1;
int N_reps = 1;


// Sprinkling Parameters
///////////////////////////////////////////
bool poisson = true;
bool make_matrix = true;
bool special = false;
bool use_transitivity = false;
bool make_sets = false;
bool make_links = false; 
const char* sets_type = "all"; 
const char* name = "cylinder";

// Shape parameters
double radius = 1;
double height = 2;
///////////////////////////////////////////

// Begin program
auto beginning = high_resolution_clock::now();

std::cout<<"\n\n============= Sprinkling into "<<name<<" ==================\n";

int iteration = 0;
for (auto dim: dims)
{
    auto start = high_resolution_clock::now();

    std::cout << "Dim = " << dim << std::endl;
    
    for (auto card : cards)
    {
        for (int rep=0; rep<N_reps; rep++)
        {
                auto repstart = high_resolution_clock::now();
                // Set up shape
                std::vector<double> center = {0.0,0.0,0.0,0.0};
                CoordinateShape shape(dim,name,center,radius,height);
                // Set up spacetime
                Spacetime S = Spacetime();
                S.BlackHoleSpacetime(dim,mass);
                // Sprinkle the causet
                SprinkledCauset C(card, S, shape, poisson,
                                make_matrix, special, use_transitivity,
                                make_sets, make_links,sets_type);

                
                std::cout << "Dim="<< dim <<", "<<(rep+1)<<"/"<<N_reps<<"\n";

                //Timing rep
                auto repend = high_resolution_clock::now();
                double duration = duration_cast<microseconds>(repend - repstart).count();
                std::cout << "Time taken generating for N = " << C._size
                << ": " << duration/pow(10,6) << " seconds" << std::endl;
        }
    
        auto mid = high_resolution_clock::now();
        double duration = duration_cast<microseconds>(mid - start).count();
        
        std::cout << "Average time taken for generating "<< N_reps
                << " causets with N = " << card << ":\n"  
                << duration/pow(10,6)/N_reps
                << " seconds\n" << std::endl;
    }   
}


auto finish = high_resolution_clock::now();
double duration = duration_cast<microseconds>(finish - beginning).count();
std::cout<<"===============================================================\n";
std::cout << "\nProgram took in total: "
        << duration/pow(10,6) << " seconds\n" << std::endl;



// end of main   
}