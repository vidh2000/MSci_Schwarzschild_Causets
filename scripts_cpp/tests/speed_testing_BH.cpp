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

#include "../causets_cpp/sprinkledcauset.h"
#include "../causets_cpp/shapes.h"
#include "../causets_cpp/spacetimes.h"

#include "../causets_cpp/functions.h"
#include "../causets_cpp/vecfunctions.h"

#include <boost/range/combine.hpp>

using namespace std::chrono;


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


// SIMULATIONS PARAMETERS (adjust only these)

int cardinality = 10000;
int dim = 4;
std::vector<int> repetitions_arr = {1};//8,7,6,5,4,3,2,

// Specify the type of causet generation
bool make_links = false; //would create future links
// Line below applies only when (make_sets||make_links)==true
const char* spacetime =  "BlackHole"; //"flat" or "BlackHole"

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////




int main(){

// Must be like that for good performance -> checked alrady
bool make_matrix = true;
const char* sets_type = "future"; //past 
bool special = false;
bool use_transitivity = false;
bool make_sets = false;
std::vector<double> center (dim, 0.0);
std::vector<int> cards = {};
std::vector<double> radii = {};
std::vector<double> durations = {};
std::vector<double> masses = {1};
for (auto mass : masses)
{
        // Make a "square" cylinder with r,h=4M, since r_S=2M
        radii.push_back(4*mass);
        durations.push_back(4*mass);
        // Keep the same density of points
        cards.push_back(cardinality*mass*mass*mass);
}

// Sprinkling Parameters
bool poisson = false;
const char* name = "cylinder";
auto beginning = high_resolution_clock::now();

std::cout<<"\n\n============ Sprinkling into "<<name<<" ===================\n";

std::cout << std::endl;


for (auto && tup : boost::combine(cards, radii, masses, durations))
{
        
        
        // Define params for causet generation
        int card;
        double radius, myduration, mass;
        boost::tie(card, radius, mass, myduration) = tup;
        std::cout << "\nSpacetime: " << spacetime << std::endl;
        std::cout << "N = " << card << ", dim = " << dim << std::endl;

        for (auto repetitions: repetitions_arr)
        {
                auto start = high_resolution_clock::now();
                // Repeat over many initialisations
                //#pragma omp parallel for
                for (int rep=0; rep<repetitions; rep++)
                {
                        //auto repstart = high_resolution_clock::now();
                        CoordinateShape shape(dim,name,center,radius,myduration);
                        Spacetime S = Spacetime();
                        if (strcmp(spacetime, "flat")==0){
                                S.FlatSpacetime(dim);}
                        else if (strcmp(spacetime, "BlackHole")==0){
                                S.BlackHoleSpacetime(dim,mass);}
                        SprinkledCauset C(card, S, shape, poisson,
                                        make_matrix, special, use_transitivity,
                                        make_sets, make_links,sets_type);

                        //Timing
                        // auto repend = high_resolution_clock::now();
                        // double duration = duration_cast<microseconds>(repend - repstart).count();
                        // std::cout << "Time taken N = " << card
                        // << ", "<<(rep+1)<<"/"<<repetitions<<": " << duration/pow(10,6)
                        // << " seconds" << std::endl;
                }
                auto mid = high_resolution_clock::now();
                double duration = duration_cast<microseconds>(mid - start).count();
                
                std::cout << "Average time taken for generating "<< repetitions <<
                        " repetitions with N = " << card << ": "
                        << duration/pow(10,6)/repetitions << " seconds" << std::endl;
        }
}

std::cout<<std::endl;
}