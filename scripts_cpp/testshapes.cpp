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


#include "causets_cpp/shapes.h"
#include "causets_cpp/functions.h"
#include "causets_cpp/vecfunctions.h"

using namespace std::chrono;
/**
 * @brief   To run this file, you must:
 *          - Go to to the folder in which it is located
 *          - Type into the command line:
 *              cd scripts_cpp
 *              g++ -g causets_cpp/shapes.cpp testshapes.cpp -std=c++17 -o testshape -O2
 *              ./testshape
 *              cd ../
 *          NOTE: running in Windows cmd does not print no matetr what
 *
 */


// Shape parameters
int card = 10000;
int dim = 4;
std::vector<double> center (dim, 0.0);


int main(){
    auto start = high_resolution_clock::now();
    //std::cout << "Starting building shape..." << std::endl;
    std::vector<const char*> names = {"ball", "bicone", "diamond", "cylinder",
                                      "cube", "cuboid"};
    for (const char* name : names)
    {
        std::cout<<"\nUSING "<<name<<std::endl;
        CoordinateShape shape(dim,name, center);
        for (auto const& p : shape._params)
            {std::cout << p.first << ' ' << p.second << '\n';}
        
        std::cout<<"Radius: "<<shape.Parameter("radius")<<std::endl;
        std::cout<<"Radius: "<<shape._params["radius"]<<std::endl;
        std::cout<<"Radius: "<<shape._params.find("radius")->second<<std::endl;
    }

    auto stop = high_resolution_clock::now();
    double duration = duration_cast<microseconds>(stop - start).count();
    std::cout << "\nTime taken: "
            << duration/pow(10,6) << " seconds" << std::endl;
    return 0;
};