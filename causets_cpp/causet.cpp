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

// KINEMATICS
double ord_fr(Causet A,
            const char* denominator, // = "choose",
            bool isdisjoined) // = true);
{

}
double ord_fr(vector<vector<int8_t>> A,
                const char* denominator,// = "choose",
                bool isdisjoined)// = true);
{

}
double ord_fr(vector<set<int>> A_future, 
                vector<set<int>> A_past,
                const char* denominator,// = "choose",
                bool isdisjoined)// = true);
{

}
double ord_fr(int a, int b,
                const char* denominator,// = "choose",
                bool isdisjoined)// = true)
{

}

static double optimiser_placeholder(){}
// typedef double (*func)();//need to pick optimiser;

template <typename F>
double MMdim_est(const char* method,// = "random",
                int d0,// = 2,
                int Nsamples,// = 20,
                int size_min,// = 10,
                double size_max,// = nan(""),
                F optimiser)// = Causet::optimiser_placeholder)
{

}   
