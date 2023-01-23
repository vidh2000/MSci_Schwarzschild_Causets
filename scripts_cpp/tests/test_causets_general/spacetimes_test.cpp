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


#include "../causets_cpp/spacetimes.h"
#include "../causets_cpp/functions.h"
#include "../causets_cpp/vecfunctions.h"

using namespace std::chrono;
/**
 * @brief   To run this file, you must:
 *          - Go to to the folder in which it is located
 *          - Type into the command line:
 *              cd scripts_cpp
 *              g++ -g ../causets_cpp/spacetimes.cpp testspacetimes.cpp -std=c++17 -o testspacetimes -O2
 *              ./testspacetimes
 *              rm testspacetimes.exe
 *              cd ../
 *          NOTE: running in Windows cmd does not print no matetr what
 *
 */


// Spacetime parameters
int card = 10000;
int dim = 4;
std::vector<double> period (dim-1, 0.0);


int main(){
    auto start = high_resolution_clock::now();
    
    std::cout<<"\n===========FLATSPACETIME==========\n";
    
    //Periodicity
    Spacetime FS = Spacetime();
    FS.FlatSpacetime(dim, period);
    std::cout<<"Should have any period? "
             <<(std::count(period.begin(), period.end(), 0.0)!=dim-1)
             <<"\nDoes have any period? "
             <<FS._isPeriodic
             <<std::endl;
    
    //Causality
    auto aretimelike = FS.Causality();
    std::vector<double> ovec (dim, 0.0);
    std::vector<double> xvec (dim, 0.0);
    std::vector<double> yvec (dim, 0.0);
    xvec[1] += 1;
    yvec[0] += 1;
    std::cout<<"Should be 0,0,0: ";
    print_vector(aretimelike(ovec, xvec, {}, 0));
    std::cout<<"Should be 1,1,0: ";
    print_vector(aretimelike(ovec, yvec, {}, 0));
    std::cout<<"Should be 1,0,1: ";
    print_vector(aretimelike(yvec, ovec, {}, 0));


    auto stop = high_resolution_clock::now();
    double duration = duration_cast<microseconds>(stop - start).count();
    std::cout << "\nTime taken: "
            << duration/pow(10,6) << " seconds" << std::endl;
    return 0;
};