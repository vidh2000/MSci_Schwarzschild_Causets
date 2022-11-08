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

#include "causets_cpp/sprinkledcauset.h"
#include "causets_cpp/shapes.h"
#include "causets_cpp/spacetimes.h"

/**
 * @brief   To run this file, you must:
 *          - Go to to the folder in which it is located
 *          - Type into the command line:
 *              g++ -g causets_cpp/causet.cpp causets_cpp/shapes.cpp
 *              causets_cpp/spacetimes.cpp causets_cpp/embeddedcauset.cpp
 *              causets_cpp/sprinkledcauset.cpp sprinkle.cpp
 *              -std=c++17 -o sprinkle -O2
 *          - This creates the executable sprinkle.exe
 *          - Run sprinkle.exe by typing into command line:
 *              ./sprinkle
 * 
 *  All individual file inside causets_cpp work.
 *  Running this script now gives:
 *   error: invalid conversion from 'CoordinateShape (*)()' to 'int'
 *                                                    [-fpermissive]
 */


// Sprinkled causet parameters
int card = 10;
int dim = 4;
FlatSpacetime S(dim);
CoordinateShape shape;
SprinkledCauset C(card,S,shape);

int main(){
    int a = 5;
    int b = 20;
    std::cout << "a x b = " << a*b << std::endl;
    std::cout << "This file works!" << std::endl;
    return 0;
};