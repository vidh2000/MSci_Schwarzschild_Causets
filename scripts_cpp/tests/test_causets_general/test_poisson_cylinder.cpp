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

using namespace std::chrono;//can't use duration


// Sprinkled causet parameters
std::vector<int> cards (100, 10000);
std::vector<const char*> names = {"cylinder"};// "ball", "cylinder", "cube", "cuboid"};
int dim = 4;
std::vector<double> center (dim, 0.0);
double hollow = 0.9;
double coeff = 4 / 3 * 3.1415926535 * (1-hollow*hollow*hollow);

double edge = 4;
std::vector<double> edges = {1,2,3,4};

// For BH
double mass = 0.5;


// Sprinkling Parameters
bool poisson = true;
bool make_matrix = false;
bool special = false;
bool use_transitivity = false;
bool make_sets = false;
bool make_links = false; 
const char* sets_type = ""; 


int main(){
    
    //std::cout << "Starting building shape..." << std::endl;
    edges.resize(dim);
    //card = (want_coords || want_matrix)? 10 : card;
    for (const char* name : names)
    {
    int counter = 0;
    for (int card : cards)
    {
        double myduration = 0.5;// std::pow(card/coeff, 0.5);
        double radius = 4;//std::pow(myduration, 0.33333333);

        auto start = high_resolution_clock::now();

        CoordinateShape shape(dim,name,center,radius,myduration,hollow,edge,edges);
        std::cout<<"\n\n============= USING "<<name<<" ====================\n";
        

        Spacetime S = Spacetime();
        S.BlackHoleSpacetime(dim);
        std::cout<<"\nSet Spacetime\n";

        SprinkledCauset C(card, S, shape, poisson,
                          make_matrix, special, use_transitivity,
                          make_sets, make_links,sets_type);
        std::cout<<"\nSprinkled\n";

        std::stringstream rstream;
        rstream << std::fixed << std::setprecision(2) << radius;
        std::string redge_s = rstream.str();
        
        std::stringstream hstream;
        hstream << std::fixed << std::setprecision(2) << hollow;
        std::string h_s = hstream.str();
        
        std::stringstream cstream;
        cstream << std::fixed << counter;
        std::string c_s = cstream.str();

        std::string filename = "../../../scripts_py/tests/test_poisson_cylinder_"
                            + std::to_string(dim)
                            + "D_N" + std::to_string(card)
                            + "_redge" + redge_s
                            + "_h" + h_s
                            + "_c" + c_s + ".txt";
        const char* path_file = filename.c_str();
        C.save_causet(path_file);
        counter += 1;


        auto stop = high_resolution_clock::now();
        double duration = duration_cast<microseconds>(stop - start).count();
        std::cout << "Time taken for N= " << card << ", dim = "<< dim << ": "
                << duration/pow(10,6) << " seconds" << std::endl;
    }
    }


    std::cout << "\n =========This file works!==========\n" << std::endl;


    
    
    return 0;
};