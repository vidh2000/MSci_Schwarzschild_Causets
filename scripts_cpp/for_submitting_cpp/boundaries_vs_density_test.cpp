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
//#include <pwd.h>

using namespace std::chrono;

std::vector<double> densities = {1,5,10,20,50,100,200,400}; //only do rho>1
int N_reps = 100;


const double pi = 3.14159265358979323846;

/**
 * @brief Add new results to previous, whereeach results is given by a map
 * <int, vector<int>>, such as for lambdas.
 * 
 * @param all_results 
 * @param newresults 
 * @param pastkeys vector<int> : do not pass an empty vector
 * @param N int : number of PREVIOUS repetitions
 */
template <typename num>
void update_distr(std::map<int, std::vector<num>> &all_results, 
                  std::map<int, num> &newresults,
                  std::vector<int> &pastkeys,
                  int N = 0)
{        
        //Create newkeys from new results
        std::vector<int> newkeys;
        for(auto it = newresults.begin(); it != newresults.end(); ++it) 
        {
            newkeys.push_back(it->first);
        }

        // Extend pastkeys if newkeys contain keys that are not in pastkeys
        // and update results with 0s in N previous rounds for new keys
        auto pastkeymax_it = std::max_element(pastkeys.begin(), pastkeys.end());
        auto newkeymax_it  = std::max_element(newkeys.begin(), newkeys.end());
        int pastkeymax = *pastkeymax_it;
        int newkeymax = *newkeymax_it;
        if (newkeymax > pastkeymax)
        {
            for (int i = pastkeymax+1; i<=newkeymax; i++)
            {
                pastkeys.push_back(i);
                all_results[i] = {};
                all_results[i].resize(N, 0);
            }
        }
        
        //Update past results with new ones
        for (int key : pastkeys)
        {
            //if newresults didn't have such key, newresults[key]=0 so no prob
            all_results[key].push_back(newresults[key]);
        }
}


/**
 * @brief get {avgs, stds}, where avgs and stds are maps int->double, where
 * int is the key (size of lambda) and double is the avg/std of the results
 * (the counts of the lambdas of that size).
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
    std::map<int, double> avgs;
    std::map<int, double> stds;

    //Get avg and std; need to avoid possible nans
    for (auto pair : all_results)
    {
        int key = pair.first;
        std::vector<num> values = pair.second;

        double sum = 0.0;
        double N = 0.0;
        std::for_each(std::begin(values), std::end(values),
                        [&](double v){if (!std::isnan(v)) {sum += v; N+=1;}}
                        );
        double avg_i = sum/N;

        double accum = 0.0;
        std::for_each(std::begin(values), std::end(values),
                        [&](double v){if (!std::isnan(v))
                        {accum += (v - avg_i) * (v - avg_i);}}
                        );
        double std_i = sqrt(accum / (N-1));

        avgs[key] = avg_i;
        stds[key] = std_i;
    }

    return {avgs, stds};
}


int main(){


int dim = 4;
double mass = 2.0;
std::vector<int> cards = {};
std::vector<double> radii = {};
std::vector<double> hollow_vals = {};
std::vector<double> durations = {};
std::vector<int> repetitions_arr = {};


std::cout << "Checking boundaries for Rho = " << std::endl;
print_vector(densities);
std::cout << "These yield causets with cardinalities:" << std::endl;

for (auto rho : densities)
{
        // Make a 4-cylinder which "just" includes all relevant links; r_S=2M
        double r_out; 
        double r_in;
        double dur;

        if (rho<=10)
        {
            dur = 3;
            r_out = 2*mass+2;
            r_in = 0;
        }   
        else if (rho>10 && rho<80)
        {
            dur = 2;
            r_out = 2*mass+1.5;
            r_in = 2*mass-1.5;
        }
        else if (rho>80)
        {
            dur = 1;
            r_out = 2*mass+1;
            r_in = 2*mass-1;
        }
        else{
            std::cout << "What the hell did you put for the density?!\n";
        }

        radii.push_back(r_out);
        hollow_vals.push_back(0);
        durations.push_back(dur);
        // Keep the same density of points, i.e such that N(M=1)=N_multiplier
        int card = rho*4./3.*pi*(r_out*r_out*r_out-r_in*r_in*r_in)*dur;
        cards.push_back(card);
        std::cout << "Cardinality = " << card << std::endl;
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

std::cout << "\nThis file counts lambdas at the constant mass and" <<
            " finds the boundary values (min/max) of radius and duration" <<
            " at different densities i.e different N_multiplier. \n \n";
std::cout << "Mass = " << mass << "\n \n";

// Variables for storage of information from each iteration
std::map<int, std::vector<double>> all_lambda_results;

int iteration = 0;
for (auto && tup : boost::combine(cards, radii, hollow_vals,
                                    densities, durations, repetitions_arr))
{
    iteration++;
    auto start = high_resolution_clock::now();

    // Store N_links for each repetition in a single iteration
    std::map<int, std::vector<double>> all_iter_lambda_results;
    std::vector<int> iter_pastkeys = {-3, -2, -1, 1};
    
    // Define params for causet generation
    int card,repetitions;
    double radius, myduration, hollow, density;
    boost::tie(card, radius, hollow, density, myduration, repetitions) = tup;
    std::cout << "======================================================\n";
    std::cout << "Main interation count: " << (iteration)<<"/"<< densities.size()
    << "\nRho = " << density << 
    ". N = " << card << ". r_S = " << 2*mass << ". Radius = " << radius <<
    ". Hollow = " << hollow <<
    ". Dimension = " << dim << ". Height = " << myduration <<
    ". Centered at t = " << -myduration/2 << ".\n";
    
    // Repeat over many causets with same parameters
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
            std::cout << "----------------------------------------------------\n";
            std::cout << "\nRho="<<density<<", "<<(rep+1)<<"/"<<repetitions<<"\n";
            std::cout << "Time taken generating for N = " << C._size
            << ": " << duration/pow(10,6) << " seconds" << std::endl;

            // Count lambdas and update results
            auto countstart = high_resolution_clock::now();
            double t_f = 0;
            std::cout << "Count lambdas starting to run now.." << std::endl;
            std::map<int, double> lambdas_distr = C.count_lambdas(t_f,2*mass);
            std::cout << "Update distr starting to run now.." << std::endl;
            update_distr(all_iter_lambda_results, lambdas_distr, 
                        iter_pastkeys, rep);

            //Timing link counting
            auto countend = high_resolution_clock::now();
            double durationlinks = duration_cast<microseconds>(
                                    countend - countstart).count();
            std::cout << "Time taken in count_lambdas for N = "
            << C._size << ": " << durationlinks/pow(10,6) << " seconds"
            << std::endl;     
    }

    std::vector<std::map<int, double>> iter_results = avg_distr(
                                                    all_iter_lambda_results, 
                                                    iter_pastkeys);
    std::map<int, double> iter_avgs = iter_results[0];
    std::map<int, double> iter_stds = iter_results[1];

    auto mid = high_resolution_clock::now();
    double duration = duration_cast<microseconds>(mid - start).count();
    
    std::cout << "Average time taken for generating "<< repetitions
            << " causets with N = " << card << ":\n"  
            << duration/pow(10,6)/repetitions
            << " seconds\n" << std::endl;   


    /* Save the average and standard deviation of lambdas for current settings
    into a text file to be read afterwards*/
    const char* homeDir = getenv("HOME");
    std::stringstream stream1;
    stream1 << std::fixed << std::setprecision(2) << mass;
    std::string mass_str = stream1.str();
    std::stringstream stream2;
    stream2 << std::fixed << std::setprecision(1) << radius;
    std::string radius_str = stream2.str();
    std::stringstream stream3;
    stream3 << std::fixed << std::setprecision(1) << myduration;
    std::string dur_str = stream3.str();
    std::stringstream stream4;
    stream4 << std::fixed << std::setprecision(0) << density;
    std::string rho_str = stream4.str();
    std::stringstream stream5;
    stream5 << std::fixed << std::setprecision(2) << hollow;
    std::string hollow_str = stream5.str();

    std::string filename = std::string(homeDir) 
                            + "/MSci_Schwarzschild_Causets/data/test_boundary_vs_density/"
                            + "Rho=" + rho_str
                            + "M=" + mass_str
                            //+ "_Nmult=" + std::to_string(N_multiplier)
                            + "_Card=" + std::to_string(card)
                            + "_r=" + radius_str
                            + "_hollow=" + hollow_str
                            + "_dur=" + dur_str
                            + ".txt";
    std::cout<<"\n==========================================================\n";
    std::cout << "Saving Iteration "<< (iteration)<<"/"<< densities.size()<<
                " to the file:\n" << filename << std::endl;  
    std::cout<<"Pastkeys found are";
    print_vector(iter_pastkeys);
    

    /////////////////////////////////////////////////////////////////////
    // CREATE/OPEN TXT FILE TO UPDATE INFORMATION

    // 1. Get the data from previous runs
    std::vector<std::string> previous_lines;
    std::string line;
    std::ifstream prev_file(filename);
    while (std::getline(prev_file, line)) 
    {
        previous_lines.push_back(line);
    }
    prev_file.close();

    // 2.1 If file didn't exist, write for first time
    if (previous_lines.size()==0)
    {
        std::ofstream out(filename);
        out<<"Nreps,"<<repetitions<<","<<std::endl;
        for (int key : iter_pastkeys)
        {
            out<<key<<"avg,"<<iter_avgs[key]<<","<<std::endl;
            out<<key<<"std,"<<iter_stds[key]<<","<<std::endl;
        }
        out.close();   
    }

    // 2.2 If file existed, rewrite file with old and new info
    else
    {
        std::ofstream out(filename);

        int N_prev_rounds = 0;
        std::string line_0 = previous_lines[0];
        for (int j = 0; j<line_0.size(); j++)
        {
            if (line_0[j] == ',')
            N_prev_rounds += 1;
        }
        N_prev_rounds -= 1; //for the , at the end of line

        for (int i = 0; i < previous_lines.size(); i++)
        {
            std::string line_i = previous_lines[i];
            if (i == 0)//Nreps line
            {
                out<<line_i<<repetitions<<","<<std::endl;
            }
            else if (i % 2) //odd i -> avg line
            {
                int key_index = i/2;
                int key = iter_pastkeys[key_index];
                out<<line_i<<iter_avgs[key]<<","<<std::endl;
            }
            else //Even i -> std line
            {
                int key_index = i/2 -1;
                int key = iter_pastkeys[key_index];
                out<<line_i<<iter_stds[key]<<","<<std::endl;
            }
        }

        //if previous file had less keys (no bigger lambda ever popped out)
        int n_prev_keys = (previous_lines.size() - 1)/2;
        if (n_prev_keys < iter_pastkeys.size())
        {
            for (int i = n_prev_keys; i<iter_pastkeys.size(); i++)
            {
                int key = iter_pastkeys[i];
                out<<key<<"avg,";
                for (int j = 0; j<N_prev_rounds; j++)
                    {out<<0<<",";}
                out<<iter_avgs[key]<<","<<std::endl;
                out<<key<<"std,";
                for (int j = 0; j<N_prev_rounds; j++)
                    {out<<0<<",";}
                out<<iter_stds[key]<<","<<std::endl;
            }
        }
        out.close();
    }
}


auto finish = high_resolution_clock::now();
double duration = duration_cast<microseconds>(finish - beginning).count();
std::cout<<"\n=============================================================\n";
std::cout << "\nProgram took in total: "
        << duration/pow(10,6) << " seconds\n" << std::endl;

}