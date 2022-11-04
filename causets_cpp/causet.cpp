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

using namespace std::chrono;


#include "functions.h"
//#include "D:\Documents\Sola\Imperial College London\Year 4\MSci project\Project\causets_code\causets_cpp\MyVecFunctions.h"
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

// CONSTRUCTOR
Causet::Causet(vector<vector<double>> coordinates,
                const char* method)
{
    // Attributes
    coords = coordinates;
    size = coords.size();
    dim = coords[0].size();

    // Sort coordinates by time - natural label
    std::sort(coords.begin(), coords.end(),
            [](auto const& lhs, auto const& rhs)
        {return lhs[0] < rhs[0];});

    /*
    Choose whether you'll build the causet as a matrix or
    specify it as a vector of sets, each set representing the element's past
    */
    if (method == "pasts"){
        make_pasts();
    }
    else if (method == "cmatrix"){
        make_cmatrix();
    }
    else{
        std::cout << "method options are: 'pasts' or 'cmatrix'" << std::endl;
    }
};

// Methods of constructing the causal set



int main(){

/*vector<vector<double>> coords = {{0.1,0.3,0.2},
                                  {1,2,3},
                                  {2,4,2},
                                  {1.5,6.3,4}};
*/
int DIM = 4;
int N = 10000;
vector<vector<double>> coords = generate_2Dvector(N,DIM,0,2);
//std::cout << "This file works" << std::endl;

auto start = high_resolution_clock::now();

Causet c(coords,"cmatrix");

auto stop = high_resolution_clock::now();
double duration = duration_cast<microseconds>(stop - start).count();
std::cout << "Time taken by function in D=" << DIM << ": "
         << duration/pow(10,6) << " seconds" << std::endl;
 



//for (int i=0; i<c.size; i++){
//print_set(c.pasts[500]);


return 0;
}

/*
// KINEMATICS
int Causet::size();
int Causet::Card();
double Causet::ord_fr(Causet A,
                const char* denominator = "choose",
                bool isdisjoined = true);
double Causet::ord_fr(vector<set<int>> A_future = {}, 
                vector<set<int>> A_past = {},
                const char* denominator = "choose",
                bool isdisjoined = true);
double Causet::ord_fr(int a, int b,
                const char* denominator = "choose",
                bool isdisjoined = true);
static double Causet::optimiser_placeholder();
//typedef double (*func)();//need to pick optimiser;
double Causet::MMdim_est(const char* method = "random",
                int d0 = 2, int Nsamples = 20,
                int size_min = 10, double size_max = nan(""),
                func optimiser = Causet::optimiser_placeholder);   


// INTERVAL
Causet Causet::Interval(int a, int b,
                bool includeBoundary = true,
                bool disjoin = false,
                const char* createmethod = "set");
int Causet::IntervalCard(int a, int b, bool includeBoundary = true);

// CAUSET REPRESENTATION & SAVING (to be added...)
vector<vector<double>> Causet::CMatrix (const char* method="causality",
                                    set<int> labels = {});
vector<vector<double>> Causet::CTMatrix (const char* method="causality",
                                    set<int> labels = {});
vector<vector<double>> Causet::LMatrix (const char* method="causality",
                                    set<int> labels = {});
void Causet::saveCasCSV(const char* filename);
void Causet::saveCasTXT(const char* filename);
void Causet::saveLasCSV(const char* filename);
void Causet::saveLasTXT(const char* filename);
*/
