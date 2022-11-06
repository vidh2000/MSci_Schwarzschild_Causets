/// \authors Vid Homsak, Stefano Veroni
/// \date 31/10/2022

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <set>
#include <stack>
#include <stdio.h>
#include <stdexcept>
#include <string>
#include <vector>
#include <chrono>
#include <unordered_set>
#include <random>

#include "functions.h"
#include "MyVecFunctions.h"
#include "causet.h"

using std::vector;
using std::set;
using std::unordered_set;


/**
 * @brief Causet class.
 *
 * 
 * @param causet: a vector of vectors of integers. //not implemented
 *  Essentially it is the causal matrix where (M)_{i,j}
 *  can take value of 1 if e_j<e_i, 0 if they aren't related and
 *  -1 if e_j<e_i and they are a link.
 * @param pasts: a vector of sets for all elements.
 * 
 * @param pastLinks: a vector of sets containing past links to each element
 */

// CONSTRUCTORS
Causet::Causet(){}
Causet::Causet(vector<vector<double>> Cmatrix, 
                bool past_links, // = false,
                bool fut_links) // = false);
{

}
Causet::Causet(vector<vector<double>> coordinates,
               const char* method) // = "pasts");
{

}

Causet::Causet(vector<vector<double>> coordinates,
                const char* method)
{

}

///////////////////////////////////////////////////////////////////////////////
// ORDERING FRACTION FUNCTIONS (OVERRIDING FOR TYPES)
///////////////////////////////////////////////////////////////////////////////
/**
 * @brief   Find the ordering fraction of an interval between a and b: 
            the ratio of actual relations over possible in such interval.
 * 
 * @param mode: string
            Use as denominator:
            - 'choose' -> |A|(|A|-1)/2, i.e. |A| choose 2 (Default).
            - 'n2'     -> (|A|^2)/2.
 *            
 * @return  Ordering fraction of Alexandrov Interval
            This is nrelations / (N choose 2)
 */

double Causet::ord_fr(Causet A,
            const char* denominator, // = "choose",
            bool isdisjoined) // = true);
{
    if (_CMatrix.size())
    {
        return Causet::ord_fr(*this._Cmatrix, denominator, isdisjoined);    
    }
    else if (_pasts.size())
    {
        return Causet::ord_fr(_futures,_pasts,denominator,isdisjoined);
    }
}

double Causet::ord_fr(vector<vector<int8_t>> A,
                const char* denominator,// = "choose",
                bool isdisjoined)// = true);
{

}
double Causet::ord_fr(vector<set<int>> A_futures, 
                vector<set<int>> A_pasts,
                const char* denominator,// = "choose",
                bool isdisjoined)// = true);
{

}
double Causet::ord_fr(int a, int b,
                const char* denominator,// = "choose",
                bool isdisjoined)// = true)
{

}

///////////////////////////////////////////////////////////////////////////////
// Dimension estimator
///////////////////////////////////////////////////////////////////////////////
/**
 * @brief Use Myrheim-Meyers dimensional estimator to compute the 
          fractal dimension (not necesseraly int).
        
 * 
 * @param method: str 
            - 'random': randomly sample.
            - 'big': take all events with no past, all with no future
                     and apply estimator ro their combinations.
 *  
 * @param d0: float 
            Initial guess for dimension.
            Default is 2.
 * 
 * @param Nsamples: int 
            Times to iterate procedure to then average on if method "random".
            Default is 20.
 * 
 * @param size_min: int\n
            Minimum size of Alexandrov Sets on which you apply estimators.
            Default is 20 elements.
 * 
 * @param size_max: int\n
            Maximum size of Alexandrov Sets on which you apply estimators.
            Default and highly recommended is np.Inf.
 * 
 * @return
        - dimension estimate: float
        - dimension std: float
 */
double Causet::MM_drelation(double d)
{
    double a = std::tgamma(d+1);
    double b = std::tgamma(d/2);
    double c = 4* std::tgamma(3*d/2);
    return a*b/c;
}

vector<double> Causet::MMdim_est(const char* method,// = "random",
                int d0,// = 2,
                int Nsamples,// = 20,
                int size_min,// = 10,
                double size_max)// = nan("")
{
    std::cout << "NOTE: MMd works only in flat spacetime" << std::endl;

    auto MM_to_solve = [](double d, double ord_fr){
        return Causet::MM_drelation(d) - ord_fr/2;
        };
    
    // Variables to be used
    int* N = &_size;
    std::vector<double> destimates;
    
    if (method == "random")
    {
        int isample = 0;
        int fails = 0;
        int successes = 0;

        while (isample < Nsamples)
        {
            if ((fails>= 1000) && (successes == 0))
            {
                std::cout << "Found 0/1000 OK Alexandrov intervals. \
                Causet portion too smol. Returning Dim<0 values.";
                return {-1,-1};
            }

            // Pick two random elements

            // Define mersenne_twister_engine Random Gen. (with random seed)
            std::random_device rd;
            int seed = rd();
            std::mt19937 gen(seed);
            std::uniform_real_distribution<> dis(0,*N);
            int e1 = (int) dis(gen), e2 =(int) dis(gen);
            int* a = nullptr, int* b = nullptr;

            if (e1 == e2){
                fails += 1;
                continue;
            }
            else if (e1 < e2){
                a = &e1;
                b = &e2;
            }
            else if (e1>e2){
                a = &e2;
                b = &e1;
            }
            else{
                fails += 1;
                continue;
            }
            
            int n = IntervalCard(*a, *b);
            if () //to be continued
            {

            }
        }
    }

}   
