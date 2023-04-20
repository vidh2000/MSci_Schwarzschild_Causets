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
 *              g++ -g ../causets_cpp/causet.cpp ../causets_cpp/shapes.cpp ../causets_cpp/spacetimes.cpp ../causets_cpp/embeddedcauset.cpp ../causets_cpp/sprinkledcauset.cpp test_matrix_sets.cpp -std=c++17 -o tms -O2
 *          - This creates the executable tms.exe
 *          - Run sprinkle.exe by typing into command line:
 *              ./tms
 * 
 *          NOTE: running in Windows cmd does not print no matetr what
 * 
 */


// Sprinkled causet parameters for Part2 - MMdim
int card = 256;
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
    cout<<"\n==== CHECK MATRIX AND SETS REPRESENTATIONS ARE EQUIVALENT ====\n";
    auto start = high_resolution_clock::now();
    cout << "\n\n ====== 2. CHECK EMPIRICAL CAUSET  =======\n";
    /* 
      IT LOOKS LIKE THIS:
          6    8    9  10
           \  /\   / /
            \/  \ //
            5    7
             \  /
              \/
              4
              |
              3
              /\
             /  \
            1    2
             \  /
              \/
              0
    */
   vector<vector<double>> coords = { { 0., 0.},
                                     { 2.,-1.},
                                     { 2., 1.},

                                     { 4., 0.},
                                     { 5., 0.},
                                     { 7.,-1.},

                                     { 9.,-1.5},
                                     { 7., 1.},
                                     { 9., 0.},

                                     { 9., 1.5},
                                     { 9., 2.5}};
    vector<vector<double>> myC = {{0, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1},
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
    vector<vector<int>> mypasts = { { },
                                    {0},
                                    {0},

                                    {0,1,2},
                                    {0,1,2,3},
                                    {0,1,2,3,4},

                                    {0,1,2,3,4,5},
                                    {0,1,2,3,4},
                                    {0,1,2,3,4,5,7},
                                    
                                    {0,1,2,3,4,7},
                                    {0,1,2,3,4,7}};
    
    vector<vector<int>> mypast_links = {{ },
                                        {0},
                                        {0},

                                        {1,2},
                                        {3},
                                        {4},

                                        {5},
                                        {4},
                                        {5, 7},
                                        
                                        {7},
                                        {7}};

    vector<vector<int>> myfuts = {  {1,2,3,4,5,6,7,8,9,10},
                                    {3,4,5,6,7,8,9,10},
                                    {3,4,5,6,7,8,9,10},

                                    {4,5,6,7,8,9,10},
                                    {5,6,7,8,9,10},
                                    {6,8},

                                    { },
                                    {8,9,10},
                                    {},
                                    
                                    {},
                                    {}};
    
    vector<vector<int>> myfut_links = { {1,2},
                                        {3},
                                        {3},

                                        {4},
                                        {5,7},
                                        {6, 8},

                                        { },
                                        {8,9,10},
                                        {},
                                        
                                        {},
                                        {}};

    Spacetime Fst = Spacetime(); Fst.FlatSpacetime(2);
    CoordinateShape S = CoordinateShape(2, "cube", {}, 1, 20);
    EmbeddedCauset Cp = EmbeddedCauset(Fst, S, coords, true, false, true,
                                       true, true, "past");
    EmbeddedCauset Cf = EmbeddedCauset(Fst, S, coords, true, false, true,
                                       true, true, "future");
    EmbeddedCauset C_all = EmbeddedCauset(Fst, S, coords, true, false, true,
                                       true, true, "all with links");
    C_all.save_causet("../../data/known_causet_from_matrixSetsTest.txt");
    //A closer look to the messy ones
    cout<<"Testing some Causalities, the following should all return 0\n";
    cout<<Fst.Flat_causal(coords[1], coords[2])[0]<<endl;
    cout<<Fst.Flat_causal(coords[5], coords[7])[0]<<endl;
    cout<<Fst.Flat_causal(coords[6], coords[7])[0]<<endl;
    cout<<Fst.Flat_causal(coords[6], coords[8])[0]<<endl;
    cout<<Fst.Flat_causal(coords[6], coords[9])[0]<<endl;
    cout<<Fst.Flat_causal(coords[6], coords[10])[0]<<endl;


    cout<<"\n ----- CMATRIX ------\n";
    cout<<"Should have the following CMatrix\n";
    print_vector(myC);
    cout<<"Our CMatrix minus One done with the past (hope all 0s)\n";
    cout<<"{";
    for (int i = 0; i<myC.size(); i++)
    {
        cout<<"{";
        for (int j = 0; j<myC.size(); j++)
            {cout<<myC[i][j]-Cp._CMatrix[i][j]<<", ";}
        cout<<"}";
    }
    cout<<"}\n";
    cout<<"Our CMatrix minus One done with the future (hope all 0s)\n";
    cout<<"{";
    for (int i = 0; i<myC.size(); i++)
    {
        cout<<"{";
        for (int j = 0; j<myC.size(); j++)
            {cout<<myC[i][j]-Cf._CMatrix[i][j]<<", ";}
        cout<<"}";
    }
    cout<<"}\n";
    cout<<"Our CMatrix minus One done with all (hope all 0s)\n";
    cout<<"{";
    for (int i = 0; i<myC.size(); i++)
    {
        cout<<"{";
        for (int j = 0; j<myC.size(); j++)
            {cout<<myC[i][j]-C_all._CMatrix[i][j]<<", ";}
        cout<<"}";
    }
    cout<<"}\n";


    cout<<"\n ----- PASTS ------\n";
    cout<<"Should have the following Past Sets\n";
    print_vector(mypasts);
    cout<<"We got with past\n";
    cout<<"{";
    for (auto past_i : Cp._pasts)
    {
        cout<<"{";
        for (int e_i : past_i)
            {cout<<e_i<<", ";}
        cout<<"}";
    }
    cout<<"}\n";
    cout<<"We got with all with links\n";
    cout<<"{";
    for (auto past_i : C_all._pasts)
    {
        cout<<"{";
        for (int e_i : past_i)
            {cout<<e_i<<", ";}
        cout<<"}";
    }
    cout<<"}\n";
    cout<<"Should have the following Past Links\n";
    print_vector(mypast_links);
    cout<<"We got with pasts\n";
    cout<<"{";
    for (auto past_i : Cp._past_links)
    {
        cout<<"{";
        for (int e_i : past_i)
            {cout<<e_i<<", ";}
        cout<<"}";
    }
    cout<<"}\n";
    cout<<"We got with all with links\n";
    cout<<"{";
    for (auto past_i : C_all._past_links)
    {
        cout<<"{";
        for (int e_i : past_i)
            {cout<<e_i<<", ";}
        cout<<"}";
    }
    cout<<"}\n";


    cout<<"\n ----- FUTURES ------\n";
    cout<<"Should have the following Future Sets\n";
    print_vector(myfuts);
    cout<<"We got from future\n";
    cout<<"{";
    for (auto fut_i : Cf._futures)
    {
        cout<<"{";
        for (int e_i : fut_i)
            {cout<<e_i<<", ";}
        cout<<"}\n";
    }
    cout<<"}\n";
    cout<<"We got from all with links\n";
    cout<<"{";
    for (auto fut_i : C_all._futures)
    {
        cout<<"{";
        for (int e_i : fut_i)
            {cout<<e_i<<", ";}
        cout<<"}\n";
    }
    cout<<"}\n";
    cout<<"Should have the following Future Links\n";
    print_vector(myfut_links);
    cout<<"We got from future\n";
    cout<<"{";
    for (auto fut_i : Cf._future_links)
    {
        cout<<"{";
        for (int e_i : fut_i)
            {cout<<e_i<<", ";}
        cout<<"}\n";
    }
    cout<<"}\n";
    cout<<"We got from all with links\n";
    cout<<"{";
    for (auto fut_i : C_all._future_links)
    {
        cout<<"{";
        for (int e_i : fut_i)
            {cout<<e_i<<", ";}
        cout<<"}\n";
    }
    cout<<"}\n";
    
    
    return 0;


    cout << "\n\n ====== 2. CHECK MMDim GIVES SAME RESULT  =======\n";
    std::vector<const char*> names = {"bicone"};
    for (const char* name : names)
    {
        CoordinateShape shape(dim,name,center,radius);
        std::cout<<"\n\n======= USING "<<name<<" ===========\n";
        

        Spacetime S = Spacetime();
        S.FlatSpacetime(dim);
        SprinkledCauset C(card, S, shape, poisson,
                          make_matrix, special, use_transitivity,
                          make_sets, make_links, sets_type);

        cout << "MM estimation with matrix (mean, std):";
        vector<double> MMd_result = C.MMdim_est("big", 20,
                                    vecmin(std::vector<int> {1000,C._size/2}),
                                    card, true);
        print_vector(MMd_result);
        MMd_result = C.MMdim_est("big", 20,
                                    vecmin(std::vector<int> {1000,C._size/2}),
                                    card, false);
        cout << "MM estimation with sets (mean, std):";
        print_vector(MMd_result);
    }


    auto stop = high_resolution_clock::now();
    double duration = duration_cast<microseconds>(stop - start).count();
    cout << "Time taken: " << duration/pow(10,6) << " seconds\n";
    return 0;
};