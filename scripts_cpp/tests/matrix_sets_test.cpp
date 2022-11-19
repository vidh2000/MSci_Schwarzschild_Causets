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

using namespace std::chrono;
using std::cout;
using std::endl;
using std::vector;
/**
 * @brief   To run this file, you must:
 *          - Go to to the folder in which it is located
 *          - Type into the command line:
 *              cd scripts_cpp
 *              g++ -g ../causets_cpp/causet.cpp ../causets_cpp/shapes.cpp ../causets_cpp/spacetimes.cpp ../causets_cpp/embeddedcauset.cpp ../causets_cpp/sprinkledcauset.cpp test_matrix_sets.cpp -std=c++17 -o tms -O2
 *          - This creates the executable tms.exe
 *          - Run sprinkle.exe by typing into command line:
 *              ./tms
 * 
 *          NOTE: running in Windows cmd does not print no matetr what
 * 
 */


// Sprinkled causet parameters
int card = 1000;
int dim = 4;
std::vector<double> center (dim, 0.0);
double radius = 4.0;


bool poisson = false;
bool make_matrix = false;
bool special = true;
bool use_transitivity = false;
bool make_sets = true;
bool make_links = false;
const char* sets_type = "all"; // "both only", "all", "pasts", "futures"

int main(){
    auto start = high_resolution_clock::now();
    std::vector<const char*> names = {"bicone"};
    for (const char* name : names)
    {
        CoordinateShape shape(dim,name,center,radius);
        std::cout<<"\n\n============= USING "<<name<<" ====================\n";
        

        Spacetime S = Spacetime();
        S.FlatSpacetime(dim);
        SprinkledCauset C(card, S, shape, poisson,
                          make_matrix, special, use_transitivity,
                          make_sets, make_links, sets_type);


        cout << "\n1. CHECK MATRIX AND SETS REPRESENTATIONS ARE EQUIVALENT\n";
        cout << "MM estimation with matrix (mean, std):";
        vector<double> MMd_result = C.MMdim_est("big", 20,
                                    vecmin(std::vector<int> {1000,C._size/2}),
                                    card, true);
        print_vector(MMd_result);
        MMd_result = C.MMdim_est("big", 20,
                                    vecmin(std::vector<int> {1000,C._size/2}),
                                    card, false);
        cout << "MM estimation with matrix (mean, std):";
        print_vector(MMd_result);
    }


    std::cout << "\n =========This file works!==========\n" << std::endl;


    auto stop = high_resolution_clock::now();
    double duration = duration_cast<microseconds>(stop - start).count();
    cout << "Time taken: "
            << duration/pow(10,6) << " seconds" << std::endl;
    return 0;
};