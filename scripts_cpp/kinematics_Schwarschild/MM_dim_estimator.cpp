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
 */
void MMdim_eqn(double d, std::vector<double> C_k_arr)
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
}


int main(){

return 0;

// end of main   
}