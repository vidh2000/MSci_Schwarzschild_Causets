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


// SIMULATIONS PARAMETERS (adjust only these)
std::vector<double> masses = {1,2,3,4};
int N_multiplier = 300;
int repetitions = 5;


int main(){

int dim = 4; //want it to be "hard coded = 4"
std::vector<double> center (dim, 0.0);
std::vector<int> cards = {};
std::vector<double> radii = {};
for (auto mass : masses)
{
        // Make a "square" cylinder with r,h=4M, since r_S=2M
        radii.push_back(4*mass);
        // Keep the same density of points
        cards.push_back(N_multiplier*mass*mass);
}

// Sprinkling Parameters
bool poisson = false;
bool make_matrix = true;
bool special = false;
bool use_transitivity = true;
bool make_sets = false;
bool make_links = true;
const char* sets_type = "future";
const char* name = "cylinder";
auto beginning = high_resolution_clock::now();

std::cout<<"\n\n============ Sprinkling into "<<name<<" ===================\n";

// Variables for storage of information from each causet realisation
std::vector<double> N_links_avgs = {};
std::vector<double> N_links_stds = {};

for (auto && tup : boost::combine(cards, radii, masses))
{
        auto start = high_resolution_clock::now();
        
        // Store N_links for each repetition
        std::vector<int> N_links_arr = {};
        // Define params for causet generation
        int card;
        double radius, mass;
        boost::tie(card, radius, mass) = tup;
        std::cout << "N = " << card << ", dim = " << dim << std::endl;
        // Repeat over many initialisations
        for (int rep=0; rep<repetitions; rep++)
        {
                std::cout << "M="<<mass<<", "<<(rep+1)<<"/"<<repetitions<<"\n";
                CoordinateShape shape(dim,name,center,radius);
                Spacetime S = Spacetime();
                S.BlackHoleSpacetime(dim,mass);
                SprinkledCauset C(card, S, shape, poisson,
                                make_matrix, special, use_transitivity,
                                make_sets, make_links,sets_type);
                // Count links and store it
                double t_f = radius;
                int N_links = C.count_links(t_f);
                N_links_arr.push_back(N_links);
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
        
        //std::cout << "M = " << mass << ", N_links = " << N_links_avg
        //        << " +- " << N_links_std  << std::endl;
        std::cout << "Average time taken for generating, N = " << card
                << ", (r,h) = " << "("<<radius<<","<< myduration
                <<"):\n" << duration/pow(10,6)/repetitions << " seconds" << std::endl;
        
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