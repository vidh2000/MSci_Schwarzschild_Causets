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

//#include "causets_cpp/sprinkledcauset.h"
#include "causets_cpp/shapes.h"
//#include "causets_cpp/spacetimes.h"

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


// Sprinkled causet parameters
int card = 10;
int dim = 2;
std::vector<double> center = {0,0};


int main(){
    //FlatSpacetime S(dim);
    //std::cout << "Starting building shape..." << std::endl;
    CoordinateShape shape(dim,"bicone", center, 2.0);
    std::cout << shape._params["radius"] << std::endl;
    std::cout << "sprinkle.cpp 'shape' radius= " << shape.Parameter("radius") << std::endl;
    //std::cout << "Coordinate Shape Created!" << std::endl;
    //SprinkledCauset C(card,S,shape);
    //std::cout << "Sprinkled Causet Created!" << std::endl;

    std::cout << "\n ====================================\n \
                    This file works!" << std::endl;
    return 0;
};