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
std::vector<int> cards = {100, 200, 300, 400, 500, 750, 1000};
std::vector<const char*> names = {"cylinder"};// "ball", "cylinder", "cube", "cuboid"};
int dim = 3;
std::vector<double> center (dim, 0.0);
double radius = 1.5;
double myduration = 5;
double hollow = 0.5;
double edge = 1.5;
std::vector<double> edges = {1,2,3,4};

// For BH
double mass = 0.5;


// Sprinkling Parameters
bool poisson = false;
bool make_matrix = false;
bool special = false;
bool use_transitivity = false;
bool make_sets = true;
bool make_links = false; 
const char* sets_type = "all with links"; 


int main(){
    
    //std::cout << "Starting building shape..." << std::endl;
    edges.resize(dim);
    //card = (want_coords || want_matrix)? 10 : card;
    for (const char* name : names)
    {
    for (int card : cards)
    {
        auto start = high_resolution_clock::now();

        CoordinateShape shape(dim,name,center,radius,myduration,hollow,edge,edges);
        std::cout<<"\n\n============= USING "<<name<<" ====================\n";
        

        Spacetime S = Spacetime();
        S.BlackHoleSpacetime(dim, mass);
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

        std::string filename = "../../../data/data_for_plotting/blackhole"
                            + std::to_string(dim)
                            + "D_N" + std::to_string(card)
                            + "_redge" + redge_s
                            + "_h" + h_s
                            +  ".txt";
        const char* path_file = filename.c_str();
        C.save_causet(path_file);


        auto stop = high_resolution_clock::now();
        double duration = duration_cast<microseconds>(stop - start).count();
        std::cout << "Time taken for N= " << card << ", dim = "<< dim << ": "
                << duration/pow(10,6) << " seconds" << std::endl;
    }
    }


    std::cout << "\n =========This file works!==========\n" << std::endl;


    
    
    return 0;
};