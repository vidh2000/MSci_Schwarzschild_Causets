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

#include "causets_cpp/sprinkledcauset.h"
#include "causets_cpp/shapes.h"
#include "causets_cpp/spacetimes.h"

#include "causets_cpp/functions.h"
#include "causets_cpp/vecfunctions.h"

using namespace std::chrono;
/**
 * @brief   To run this file, you must:
 *          - Go to to the folder in which it is located
 *          - Type into the command line:
 *              cd scripts_cpp
 *              g++ -g causets_cpp/causet.cpp causets_cpp/shapes.cpp causets_cpp/spacetimes.cpp causets_cpp/embeddedcauset.cpp causets_cpp/sprinkledcauset.cpp sprinkle.cpp -std=c++17 -o sprinkle -O2
 *          - This creates the executable sprinkle.exe
 *          - Run sprinkle.exe by typing into command line:
 *              ./sprinkle
 * 
 *          NOTE: running in Windows cmd does not print no matetr what
 * 
 */
// gives 0 but doesn't create a new entry std::cout << "sprinkle.cpp 'shape' radius= " << shape._params.find("radius")->second << std::endl;
// gives 0 but creates a new entry std::cout << "sprinkle.cpp 'shape' radius= " << shape._params["radius"]<< std::endl;
   

// Sprinkled causet parameters
int card = 100;
int dim = 3;
std::vector<double> center (dim, 0.0);

bool poisson = false;
bool make_matrix = true;
bool special = true;
bool use_transitivity = true;
bool make_sets = false;
bool make_links = false;

int main(){
    auto start = high_resolution_clock::now();
    //std::cout << "Starting building shape..." << std::endl;
    std::vector<const char*> names = {"ball", "bicone", "diamond", "cylinder",
                                      "cube", "cuboid"};
    
    for (const char* name : names)
    {
        std::cout<<"\n\n============= USING "<<name<<" ====================\n";
        std::cout << "What are the just-Shape's parameters?\n";
        CoordinateShape shape(dim,name, center);
        for (auto const& p : shape._params)
            {std::cout << "-- " << p.first << "=" << p.second << '\n';}
        

        FlatSpacetime S(dim);
        SprinkledCauset C(card,S,shape,poisson,make_matrix,special,
                use_transitivity,make_sets,make_links);
        

        std::cout << "\nWhat are the Causet's Shape's parameters at the end?\n";
        for (auto const& p : C._shape._params){
        std::cout << "-- " << p.first << '=' << p.second << '\n';}
        std::cout<<"-- cardinality=" << C._size << "\n";
        std::cout<<"-- dim=" << C.spacetime_dim() << "\n";

        std::cout << "\nLet's look at the sprinkled values\n";
        std::cout<<"Max Eu Distance: "<<C.max_eu_dist()<<std::endl;
        std::cout<<"Max Sp Radius  : "<<C.max_sp_rad() <<std::endl;
        std::cout<<"Max Time       : "<<C.max_along(0) <<std::endl;
        std::cout<<"Min Time       : "<<C.min_along(0) <<std::endl;
        std::cout<<"Max Along x    : "<<C.max_along(1) <<std::endl;
        std::cout<<"Min Along x    : "<<C.min_along(1) <<std::endl;
    }

    std::cout << "\n ====================================\n \
                    This file works!" << std::endl;


    auto stop = high_resolution_clock::now();
    double duration = duration_cast<microseconds>(stop - start).count();
    std::cout << "\nTime taken: "
            << duration/pow(10,6) << " seconds" << std::endl;
    return 0;
};