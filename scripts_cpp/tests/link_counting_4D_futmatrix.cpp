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
#include <omp.h>

using namespace std::chrono;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


// SIMULATIONS PARAMETERS (adjust only these)
std::vector<double> masses = {1,1.5,2,2.5,3,3.5,4};
int N_multiplier = 400;
//int N_reps = 20;
std::vector<int> repetitions_arr = {100,100,100,100,50,20,20};

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

int main(){

int dim = 4; //want it to be "hard coded = 4"
std::vector<int> cards = {};
std::vector<double> radii = {};
std::vector<double> durations = {};
//std::vector<int> repetitions_arr = {};

for (auto mass : masses)
{
        // Make a "square" cylinder with r,h=4M, since r_S=2M
        radii.push_back(4*mass);
        durations.push_back(4*mass);
        // Keep the same density of points
        cards.push_back(N_multiplier*mass*mass*mass*mass);
        // Add # of repetitions for each mass
        //repetitions_arr.push_back(N_reps);
}

// Sprinkling Parameters
bool poisson = false;
bool make_matrix = true;
bool special = false;
bool use_transitivity = false;
bool make_sets = false;
bool make_links = false; 
const char* sets_type = "future"; 
const char* name = "cylinder";
auto beginning = high_resolution_clock::now();

std::cout<<"\n\n============ Sprinkling into "<<name<<" ===================\n";
std::cout << "Doing CMatrix and inferring future links from it\n \n";
std::cout << "N_multiplier = " << N_multiplier << "\n \n";

// Variables for storage of information from each causet realisation
std::vector<double> N_links_avgs = {};
std::vector<double> N_links_stds = {};

int iteration = 0;
for (auto && tup : boost::combine(cards, radii, masses, durations, repetitions_arr))
{
        iteration++;
        auto start = high_resolution_clock::now();

        // Store N_links for each repetition
        std::vector<int> N_links_arr = {};
        // Define params for causet generation
        int card,repetitions;
        double radius, myduration, mass;
        boost::tie(card, radius, mass, myduration, repetitions) = tup;
        std::cout << "======================================================\n";
        std::cout << "Main interation count: " << (iteration)<<"/"<< masses.size() <<
        "\nN = " << card << ". r_S = " << 2*mass << ". Radius = " << radius <<
        ". Dimension = " << dim << ". Height = " << myduration <<
        ". Centered at t = " << -myduration/2 << ".\n";
        
        // Repeat over many initialisations
        for (int rep=0; rep<repetitions; rep++)
        {
                auto repstart = high_resolution_clock::now();
                // Set up shape
                std::vector<double> center = {-myduration/2,0.0,0.0,0.0};
                CoordinateShape shape(dim,name,center,radius,myduration);
                // Set up spacetime
                Spacetime S = Spacetime();
                S.BlackHoleSpacetime(dim,mass);
                // Sprinkle the causet
                SprinkledCauset C(card, S, shape, poisson,
                                make_matrix, special, use_transitivity,
                                make_sets, make_links,sets_type);

                //Timing generation
                auto repend = high_resolution_clock::now();
                double duration = duration_cast<microseconds>(repend - repstart).count();
                std::cout << "M="<<mass<<", "<<(rep+1)<<"/"<<repetitions<<"\n";
                std::cout << "Time taken generating for N = " << card
                << ": " << duration/pow(10,6) << " seconds" << std::endl;

                // Count links and store it
                auto linkcountstart = high_resolution_clock::now();
                double t_f = 0;
                double N_links = C.count_links_fromCMatrix(t_f,2*mass)*1.0;
                N_links_arr.push_back(N_links);

                //Timing link counting
                auto linkcountend = high_resolution_clock::now();
                double durationlinks = duration_cast<microseconds>(
                        linkcountend - linkcountstart).count();
                std::cout << "N_links counted = " << N_links << std::endl;
                std::cout << "Time taken in count_links_fromCMatrix for N = "
                << card << ": " << durationlinks/pow(10,6) << " seconds"
                << std::endl;

                
        }
        double N_links_avg = mymean(N_links_arr);
        double accum = 0.0;
        std::for_each(std::begin(N_links_arr), std::end(N_links_arr),
                    [&](const double N)
                    {accum += (N - N_links_avg) * (N - N_links_avg);});

        double N_links_std = sqrt(accum / (N_links_arr.size()));
        N_links_avgs.push_back(N_links_avg);
        N_links_stds.push_back(N_links_std);
        auto mid = high_resolution_clock::now();
        double duration = duration_cast<microseconds>(mid - start).count();
        
        std::cout <<"Number of links over the horizon = " << N_links_avg
                << " +- " << N_links_std  << std::endl;
        std::cout << "Average time taken for generating "<< repetitions
                << " causets with N = " << card << ":\n"  
                << duration/pow(10,6)/repetitions
                << " seconds\n" << std::endl;
        
}

auto finish = high_resolution_clock::now();
double duration = duration_cast<microseconds>(finish - beginning).count();
std::cout << "\nProgram took in total: "
        << duration/pow(10,6) << " seconds\n" << std::endl;

std::cout << "We want N_links proportional to the mass.\nResults:" << std::endl;
for (auto && tup : boost::combine(masses, N_links_avgs, N_links_stds))
{
        double mass, N_links_avg, N_links_std;
        boost::tie(mass, N_links_avg, N_links_std) = tup;
        std::cout << "M = "<<mass<<", N_links = " << N_links_avg
                << " +- " << N_links_std  << std::endl;
}
std::cout<<std::endl;
}