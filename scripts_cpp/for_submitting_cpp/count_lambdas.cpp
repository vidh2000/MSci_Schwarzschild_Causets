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

// $HOME var get
#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>

using namespace std::chrono;

/**
 * @brief Add new results to previous, whereeach results is given by a map
 * <int, vector<int>>, such as for lambdas.
 * 
 * @param all_results 
 * @param newresults 
 * @param pastkeys 
 */
template <typename num>
void update_distr(std::map<int, std::vector<num>> &all_results, 
                  std::map<int, num> &newresults,
                  std::vector<int> &pastkeys)
{
        int N = (all_results.begin()->second).size(); //number of done reps

        //Update pastkeys
        if (pastkeys.size()==0)
        {
            for(std::map<int,std::vector<int>>::iterator it = all_results.begin();
                it != all_results.end(); ++it) 
                {
                    pastkeys.push_back(it->first);
                }
        }
        
        //Create newkeys
        std::vector<int> newkeys;
        for(std::map<int,int>::iterator it = newresults.begin(); 
            it != newresults.end(); ++it) 
        {
            newkeys.push_back(it->first);
        }

        //Update past results with new ones

        for (int key : pastkeys)
        {
            if (std::find(newkeys.begin(), newkeys.end(), key) != newkeys.end())
            {
                all_results[key].push_back(newresults[key]);
            }
            else
            {
                all_results[key].push_back(0);
            }
        }

        for (int key : newkeys)
        {
            if (std::find(pastkeys.begin(), pastkeys.end(), key) != pastkeys.end())
            {
                continue; //the case in which both contained key was treated
            }
            else
            {
                all_results[key] = {};
                all_results[key].resize(N, 0);
                all_results[key].push_back(newresults[key]);
            }
        }
}


/**
 * @brief get {avgs, stds}, where avgs and stds are maps int->double, where
 * int is the key (size of lambda) and double is the avg/std of the results (the 
 * counts of the lambdas of that size).
 * 
 * @param all_results 
 * @param newresults 
 * @param pastkeys 
 */
template <typename num>
std::vector<std::map<int, double>> avg_distr(
                                std::map<int, std::vector<num>> &all_results, 
                                std::vector<int> &pastkeys)
{
        int N = (all_results.begin()->second).size(); //number of done reps
        std::map<int, double> avgs;
        std::map<int, double> stds;

        //Get avg
        for (auto pair : all_results)
        {
            int key = pair->first;
            std::vector<num> values = pair->second;
            double avg_i = mymean(values);

            double accum = 0.0;
            std::for_each(std::begin(values), std::end(values),
                          [&](double N){accum += (N - avg_i) * (N - avg_i);}
                          );
            double std_i = sqrt(accum / N);

            avgs[key] = avg_i;
            stds[key] = std_i;
        }

        //Get std
        for (auto pair : all_results)
        {
            int key = pair->first;
            std::vector<num> values = pair->second;
            avgs[key] = mymean(values);
        }     
}



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//        PARAMETERS are inputted via COMMAND LINE in BASH script
//              - param 1 = mass (double)
//              - param 2 = N_multiplier (int)
//              - param 3 = N_reps (int)
//___________________________________________________________________________//
//////////////////////// N = N_multiplier*mass^3 //////////////////////////////
//---------------------------------------------------------------------------//
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]){

double mass = std::atof(argv[1]); 
int N_multiplier = std::atoi(argv[2]); //1000;
int N_reps = std::atoi(argv[3]);

std::cout << "PARAMETERS used in the causet generation:\n";
std::cout << "mass="<<mass<<", N_multiplier="<<N_multiplier<<", N_reps="
                <<N_reps<<"\n\n";
                
int dim = 4; //want it to be "hard coded = 4"
std::vector<double> masses = {mass};
std::vector<int> cards = {};
std::vector<double> radii = {};
std::vector<double> hollow_vals = {};
std::vector<double> durations = {};
std::vector<int> repetitions_arr = {};

for (auto mass : masses)
{
        // Make a "square" cylinder with r,h=4M, since r_S=2M
        radii.push_back(2*mass+1);
        hollow_vals.push_back((2*mass-1)/(2*mass+1));
        durations.push_back(1); // since min(t_min) ~ -3.5, 4 is adequate
        // Keep the same density of points, i.e such that N(M=1)=N_multiplier
        cards.push_back(N_multiplier/26.0*
        ((2*mass+1)*(2*mass+1)*(2*mass+1)-(2*mass-1)*(2*mass-1)*(2*mass-1)));
        // Add # of repetitions for each mass
        repetitions_arr.push_back(N_reps);
}

// Sprinkling Parameters
bool poisson = true;
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
for (auto && tup : boost::combine(cards, radii, hollow_vals,
                                         masses, durations, repetitions_arr))
{
        iteration++;
        auto start = high_resolution_clock::now();

        // Store N_links for each repetition
        std::map<int, std::vector<double>> all_lambdas_results;
        std::vector<int> pastkeys = {};

        // Define params for causet generation
        int card,repetitions;
        double radius, myduration, mass, hollow;
        boost::tie(card, radius, hollow, mass, myduration, repetitions) = tup;
        std::cout << "======================================================\n";
        std::cout << "Main interation count: " << (iteration)<<"/"<< masses.size() <<
        "\nN = " << card << ". r_S = " << 2*mass << ". Radius = " << radius <<
        ". Hollow = " << hollow <<
        ". Dimension = " << dim << ". Height = " << myduration <<
        ". Centered at t = " << -myduration/2 << ".\n";
        
        // Repeat over many initialisations
        for (int rep=0; rep<repetitions; rep++)
        {
                auto repstart = high_resolution_clock::now();
                // Set up shape
                std::vector<double> center = {-myduration/2,0.0,0.0,0.0};
                CoordinateShape shape(dim,name,center,radius,myduration,hollow);
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
                std::cout << "Time taken generating for N = " << C._size
                << ": " << duration/pow(10,6) << " seconds" << std::endl;

                // Count lambdas and store it
                auto linkcountstart = high_resolution_clock::now();
                double t_f = 0;
                std::map<int, double> lambdas_distr = C.count_lambdas(t_f,2*mass);
                update_distr(all_lambdas_results, lambdas_distr, pastkeys);

                //Timing link counting
                auto linkcountend = high_resolution_clock::now();
                double durationlinks = duration_cast<microseconds>(
                        linkcountend - linkcountstart).count();
                std::cout << "Time taken in count_lambdas for N = "
                << C._size << ": " << durationlinks/pow(10,6) << " seconds"
                << std::endl;     
        }

        std::vector<std::map<int, double>> results = avg_distr(
                                                     all_lambdas_results, 
                                                     pastkeys);
        std::map<int, double> avgs = results[0];
        std::map<int, double> stds = results[1];

        auto mid = high_resolution_clock::now();
        double duration = duration_cast<microseconds>(mid - start).count();
        
        std::cout << "Average time taken for generating "<< repetitions
                << " causets with N = " << card << ":\n"  
                << duration/pow(10,6)/repetitions
                << " seconds\n" << std::endl;     
}


auto finish = high_resolution_clock::now();
double duration = duration_cast<microseconds>(finish - beginning).count();
std::cout<<"\n=============================================================\n";
std::cout<<"===============================================================\n";
std::cout << "\nProgram took in total: "
        << duration/pow(10,6) << " seconds\n" << std::endl;

std::cout << "Parameters used:\n";
std::cout << "Dim = "<< dim << ", N_multiplier = "<< N_multiplier << std::endl;

std::cout << "Results:" << std::endl;
for (auto && tup : boost::combine(masses, N_links_avgs, N_links_stds))
{
        double mass, N_links_avg, N_links_std;
        boost::tie(mass, N_links_avg, N_links_std) = tup;
        std::cout << "M = "<<mass<< ", Card = " << cards[0]
        << ", N_links = " << N_links_avg << " +- " << N_links_std 
        << std::endl;

        /* Save the average and standard deviation of links for current settings
        into a text file to be read afterwards*/
        const char* homeDir = getenv("HOME");
        std::stringstream stream1;
        stream1 << std::fixed << std::setprecision(2) << mass;
        std::string mass_str = stream1.str();
        std::stringstream stream2;
        stream2 << std::fixed << std::setprecision(1) << radii[0];
        std::string radius_str = stream2.str();
        std::stringstream stream3;
        stream3 << std::fixed << std::setprecision(1) << durations[0];
        std::string dur_str = stream3.str();

        std::string filename = std::string(homeDir) 
                + "/MSci_Schwarzschild_Causets/data/lambdas/"
                + "M=" + mass_str
                + "_Rho=" + std::to_string(N_multiplier)
                + "_Card=" + std::to_string(cards[0])
                + "_r=" + radius_str
                + "_dur=" + dur_str
                + ".txt";
        
        std::cout << "Saving to the file: " << filename << std::endl;
        
        // Create/open the text file then write into it
        std::ofstream out(filename,std::ios_base::app);
        out << "\nN_reps, N_links_avg, N_links_std,       " <<
        N_reps << ", " << N_links_avg << ", " << N_links_std;
        out.close();

}





}