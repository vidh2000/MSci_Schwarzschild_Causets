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



// Defining constants (functions)
///////////////////////////////////////////////////////////////////////////////


/**
 * @brief Calculate n choose k and return the value as double
 * 
 * @param n Integer
 * @param k Integer
 * @return double 
 */
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
 * @brief 
 * 
 * @param d - Manifold dimension 
 * @param C_k_arr - List of C_ks (Numbers of k-long chains)
 *               for chain lengths k=1,2,3,4.
 * @return double
 */
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

int main(){




// Create causet in Schwarzschild spacetime
///////////////////////////////////////////////////////////////////////////////
                
std::vector<int> dims = {4}; 
std::vector<int> cards = {1000};
double mass = 1;
int N_reps = 10;


// Sprinkling Parameters
///////////////////////////////////////////
bool poisson = true;
bool make_matrix = true;
bool special = false;
bool use_transitivity = false;
bool make_sets = true;
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
        std::vector<double> dim_ests = {}; 

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


            // Create interval          

            // Estimate dimension of the causet
            //double d_i = estimate_MMd(C_k_arr);
            //dim_ests.push_back(d_i);

            
            std::cout << "Dim="<< dim <<", "<<(rep+1)<<"/"<<N_reps<<"\n";

            //Timing rep
            auto repend = high_resolution_clock::now();
            double duration = duration_cast<microseconds>(repend - repstart).count();
            std::cout << "Time taken generating for N = " << C._size
            << ": " << duration/pow(10,6) << " seconds" << std::endl;
        }
    
        auto mid = high_resolution_clock::now();
        double duration = duration_cast<microseconds>(mid - start).count();
        
        std::cout << "-------------------------------------\n";
        std::cout << "<t> for "<< N_reps
                << " causets in D=" << dim << " with N = " << card << ": "  
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