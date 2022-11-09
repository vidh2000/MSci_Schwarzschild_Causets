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
 *          cd scripts_cpp && g++ -g causets_cpp/causet.cpp causets_cpp/shapes.cpp causets_cpp/spacetimes.cpp causets_cpp/embeddedcauset.cpp causets_cpp/sprinkledcauset.cpp sprinkle.cpp -std=c++17 -o sprinkle -O2 && start sprinkle.exe & cd ../

 * 
 */
// gives 0 but doesn't create a new entry std::cout << "sprinkle.cpp 'shape' radius= " << shape._params.find("radius")->second << std::endl;
// gives 0 but creates a new entry std::cout << "sprinkle.cpp 'shape' radius= " << shape._params["radius"]<< std::endl;
   

// Sprinkled causet parameters
int card = 100;
int dim = 4;
std::vector<double> center(dim,0.0);

bool poisson = false;
bool make_matrix = true;
bool special = true;
bool use_transitivity = true;
bool make_sets = false;
bool make_links = false;

int main(){
    auto start = high_resolution_clock::now();
    //std::cout << "Starting building shape..." << std::endl;
    CoordinateShape shape(dim,"ball",center,2.0);

    FlatSpacetime S(dim);
    SprinkledCauset C(card,S,shape,poisson,make_matrix,special,
                use_transitivity,make_sets,make_links);
    std::cout << "\nCoordinates:\n";
    print_vector(C._coords);
    
    std::cout << "Causet Parameters:\n";
    for (auto const& p : C._shape._params){
    std::cout << "-- " << p.first << '=' << p.second << '\n';}
    std::cout<<"-- cardinality=" << C._size << "\n";
    std::cout<<"-- dim=" << C.spacetime_dim() << "\n";


    //std::cout<<"Causal Matrix:\n";
    //print_vector(C._CMatrix);
    
    std::cout<<"\n=========================================================\n"
                  <<  "This file works!" << std::endl;

    auto stop = high_resolution_clock::now();
    double duration = duration_cast<microseconds>(stop - start).count();
    std::cout << "Time taken: "
            << duration/pow(10,6) << " seconds" << std::endl;

    return 0;
};