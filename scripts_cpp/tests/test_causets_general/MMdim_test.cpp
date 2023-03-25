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

#include "../../causets_cpp/sprinkledcauset.h"
#include "../../causets_cpp/shapes.h"
#include "../../causets_cpp/spacetimes.h"
#include "../../causets_cpp/sprinkledcauset.h"
#include "../../causets_cpp/shapes.h"
#include "../../causets_cpp/spacetimes.h"

#include "../../causets_cpp/functions.h"
#include "../../causets_cpp/vecfunctions.h"
#include "../../causets_cpp/functions.h"
#include "../../causets_cpp/vecfunctions.h"

using namespace std::chrono;
using std::cout;
using std::endl;
using std::vector;

/**
 * @brief   To run this file, you must:
 *          - Go to to the folder in which it is located
 *          - Type into the command line:
 *              cd scripts_cpp
 *              g++ -g ../causets_cpp/causet.cpp ../causets_cpp/shapes.cpp ../causets_cpp/spacetimes.cpp ../causets_cpp/embeddedcauset.cpp ../causets_cpp/sprinkledcauset.cpp test_MMdim.cpp -std=c++17 -o tMM -O2
 *          - This creates the executable tMM.exe
 *          - Run sprinkle.exe by typing into command line:
 *              ./tMM
 * 
 *          NOTE: running in Windows cmd does not print no matetr what
 * 
 */


// Sprinkled causet parameters
int card = 2000;
int dim = 4;
std::vector<double> center (dim, 0.0);
double radius = 200.0;
double myduration = 5;


bool poisson = true;
bool make_matrix = true;
bool special = false;
bool use_transitivity = false;
bool make_sets = false;
bool make_links = false;
const char* sets_type = "all"; // "both only", "all", "pasts", "futures"

int main(){
    auto start = high_resolution_clock::now();
    
    cout<<"\n==============================================\n";
    cout<<"========FIRST TEST ORDERING FRACTION==========\n";
    cout<<"==============================================\n";

    cout<<"\nTEST1.1: 1D CHAIN CAUSET\n";
    int N = 8;
    vector<vector<int>> Cm (N, vector<int> (8));
    for (int i = 0; i<N-1; i++)
    {
        for (int j = i+1; j<N; j++) Cm[i][j] = 1;
    }
    Causet C = Causet(Cm);
    C.make_sets_fromC();
    cout<<"A1 has ord_fr from CMatrix = "<<C.ord_fr(3,6)<<" = 1?\n";
    cout<<"A1 has ord_fr from sets    = "<<C.ord_fr(3,6,"choose",false)<<" = 1?\n";
    cout<<"A1 has ord_fr from CMatrix = "<<C.ord_fr(0,N-1)<<" = 1?\n";
    cout<<"A1 has ord_fr from sets    = "<<C.ord_fr(3,6,"choose",false)<<" = 1?\n";


    cout<<"\nTEST1.2: A GENERAL CAUSET\n";
    /* 
      IT LOOKS LIKE THIS:
          6   8  9  10
          \  / \ | /
            5    7
             \  /
               4
               |
               3
             /   \
            1     2
             \   /
               0
    */
    Cm = {{0, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1},
          {0, 0, 0,  1, 1, 1,  1, 1, 1,  1, 1},
          {0, 0, 0,  1, 1, 1,  1, 1, 1,  1, 1},

          {0, 0, 0,  0, 1, 1,  1, 1, 1,  1, 1},
          {0, 0, 0,  0, 0, 1,  1, 1, 1,  1, 1},
          {0, 0, 0,  0, 0, 0,  1, 0, 1,  0, 0},

          {0, 0, 0,  0, 0, 0,  0, 0, 0,  0, 0},
          {0, 0, 0,  0, 0, 0,  0, 0, 1,  1, 1},
          {0, 0, 0,  0, 0, 0,  0, 0, 0,  0, 0},

          {0, 0, 0,  0, 0, 0,  0, 0, 0,  0, 0},
          {0, 0, 0,  0, 0, 0,  0, 0, 0,  0, 0}};
    C = Causet(Cm);
    cout<<"Has ord_fr from CMatrix = "<< C.ord_fr(3,8)<<endl;
    cout<<"Has ord_fr from sets    = "<< C.ord_fr(3,8,"choose",true)<<endl;
    cout<<"Should have 0.9\n";
    cout<<"Has ord_fr from CMatrix = "<< C.ord_fr(2,8)<<endl;
    cout<<"Has ord_fr from sets    = "<< C.ord_fr(2,8,"choose",true)<<endl;
    cout<<"Should have 0.93\n";
    cout<<"Has ord_fr from CMatrix = "<< C.ord_fr(3,9)<<endl;
    cout<<"Has ord_fr from sets    = "<< C.ord_fr(3,9,"choose",true)<<endl;
    cout<<"Should have 1\n";
    
    cout<<"\n========================================================\n";
    cout<<"==== SECOND TEST MYRHEIM MEYER DIMENSIONAL ESTIMATOR ====";
    cout<<"\n========================================================\n";
    Spacetime S = Spacetime();
    S.BlackHoleSpacetime(dim);
    CoordinateShape shape(dim,"cylinder",center,radius,myduration);
    // S.FlatSpacetime();
    // CoordinateShape shape(dim,"bicone",center,radius);


    cout << "\n1. CHECK MMDIM\n";
    for (int i = 0; i<10; i++){
    SprinkledCauset Cs(card, S, shape, poisson,
                        make_matrix, special, use_transitivity,
                        make_sets, make_links, sets_type);
    vector<double> MMd_result = Cs.MMdim_est("big", 10, card/100, card, true);
    cout << "\nMM estimation with CMatrix (mean, std):";
    print_vector(MMd_result);
    MMd_result = Cs.MMdim_est("big", 10, card/100, card, false);
    cout << "MM estimation with SETS   (mean, std):";
    MMd_result = Cs.MMdim_est("random", 10, card/100, card, true);
    cout << "\nMM estimation with CMatrix (mean, std):";
    print_vector(MMd_result);
    MMd_result = Cs.MMdim_est("random", 10, card/100, card, false);
    cout << "MM estimation with SETS   (mean, std):";
    print_vector(MMd_result);
    }


    auto stop = high_resolution_clock::now();
    double duration = duration_cast<microseconds>(stop - start).count();
    cout << "\n=============Time taken: "
            << duration/pow(10,6) << " seconds===============\n" << std::endl;
    return 0;
};