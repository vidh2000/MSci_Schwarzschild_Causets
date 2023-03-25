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

///////////////////////////////////////////////////////////////////////////////
// USEFUL FUNCTIONS
///////////////////////////////////////////////////////////////////////////////

template <typename num>
void update_distr(std::map<int, std::vector<num>> &all_results, 
                  std::map<int, num> &newresults,
                  std::vector<int> &pastkeys,
                  int N = 0);

template <typename num>
std::vector<std::map<int, double>> avg_distr(
                                std::map<int, std::vector<num>> &all_results, 
                                std::vector<int> &pastkeys);



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//        PARAMETERS are inputted via COMMAND LINE in BASH script
//              - param 1 = mass (double)
//              - param 2 = Rho (int)
//              - param 3 = N_reps (int)
//___________________________________________________________________________//
//////////////////////// N = Rho*mass^3 //////////////////////////////
//---------------------------------------------------------------------------//
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]){

double mass = std::atof(argv[1]); 
double Rho = std::atoi(argv[2]); //1000;
int N_reps = std::atoi(argv[3]);


std::cout << "PARAMETERS used in the causet generation:\n";
std::cout << "mass="<<mass<<", Rho="<<Rho<<", N_reps="
                <<N_reps<<"\n\n";
                
int dim = 4; //want it to be "hard coded = 4"
std::vector<double> masses = {mass};
std::vector<int> cards = {};
std::vector<double> radii = {};
std::vector<double> hollow_vals = {};
std::vector<double> durations = {};
std::vector<int> repetitions_arr = {};


// Shape Parameters
double scale = std::pow(Rho, -1.0/4.0);
std::cout << "Scale = " << scale << std::endl;
double R = 2*mass+3*scale;
std::cout << "R done"<< std::endl;
double r = 2*mass-3*scale;
std::cout << "r done"<< std::endl;
double T = 4*scale;
std::cout << "T done"<< std::endl;
double h = r/R;
std::cout << "h done"<< std::endl;
int N = Rho * (4*3.1415/3) * (R*R*R-r*r*r) * T;
std::cout << "N done"<< std::endl;
radii.push_back(R);
std::cout << "r pushed back"<< std::endl;
hollow_vals.push_back(h);
std::cout << "h pushed back"<< std::endl;
durations.push_back(T); 
std::cout << "T pushed back"<< std::endl;
cards.push_back(N);
std::cout << "N pushed back"<< std::endl;
// Add # of repetitions for each mass
repetitions_arr.push_back(N_reps);
std::cout << "N_reps pushed back"<< std::endl;


// Sprinkling Parameters
bool poisson = true;
bool make_matrix = true;
bool special = false;
bool use_transitivity = false;
bool make_sets = false;
bool make_links = false; 
const char* sets_type = "future"; 
const char* name = "cylinder";
std::cout << "Sprinkling Parameters Set"<< std::endl;
auto beginning = high_resolution_clock::now();
std::cout << "Clock Begun"<< std::endl;

std::cout<<"\n\n============ Sprinkling into "<<name<<" ===================\n";
std::cout << "Doing CMatrix and inferring future links from it\n \n";
std::cout << "Rho = " << Rho << "\n \n";

// Variables for storage of information from each iteration
std::map<int, std::vector<double>> all_lambda_results;

int iteration = 0;
for (auto && tup : boost::combine(cards, radii, hollow_vals,
                                    masses, durations, repetitions_arr))//zip
{
    iteration++;
    auto start = high_resolution_clock::now();

    // Store N_links for each repetition in a single iteration
    std::map<int, std::vector<double>> all_iter_lambda_results;
    std::vector<int> iter_pastkeys = {-3, -2, -1, 1};
    
    // Define params for causet generation
    int card,repetitions;
    double radius, myduration, mass, hollow;
    boost::tie(card, radius, hollow, mass, myduration, repetitions) = tup;
    std::cout << "======================================================\n";
    std::cout << "Main interation count: " << (iteration)<<"/"<< masses.size()
        <<"\nN = " << card << ". r_S = " << 2*mass << ". Radius = " << radius <<
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
            std::cout << "M="<<mass<<", "<<(rep+1)<<"/"<<repetitions<<"\n";
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

            // std::cout<< "OVERALL MEMORY CONSUMPTION\n";
            // std::cout<< (sizeof(C) + sizeof(lambdas_distr) 
            //           + sizeof(S) + sizeof(shape)
            //           + sizeof(all_iter_lambda_results)
            //           + sizeof(iter_pastkeys))*8
            //           <<" bits"
            //           <<std::endl;   
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
    std::stringstream stream0;
    stream0 << std::fixed << std::setprecision(2) << Rho;
    std::string rho_str = stream0.str();
    std::stringstream stream1;
    stream1 << std::fixed << std::setprecision(2) << mass;
    std::string mass_str = stream1.str();
    std::stringstream stream2;
    stream2 << std::fixed << std::setprecision(2) << radii[0];
    std::string radius_str = stream2.str();
    std::stringstream stream3;
    stream3 << std::fixed << std::setprecision(2) << durations[0];
    std::string dur_str = stream3.str();
    std::stringstream stream4;
    stream4 << std::fixed << std::setprecision(2) << hollow;
    std::string hollow_str = stream4.str();

    std::string filename = std::string(homeDir) 
                            + "/MSci_Schwarzschild_Causets/data/lambdas/"
                            + "M=" + mass_str
                            + "_Rho=" + rho_str 
                            + "_Card=" + std::to_string(cards[0])
                            + "_r=" + radius_str
                            + "_hollow=" + hollow_str
                            + "_dur=" + dur_str
                            + ".txt";
    std::cout<<"\n==========================================================\n";
    std::cout << "Saving Iteration "<< iteration <<
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
        // First line: labels
        out<<std::left<<std::setw(11);
        out<<"Nreps,";
        for (int key : iter_pastkeys)
        {
            out<<std::right<<std::setw(3);
            out<<key;
            out<<std::left<<std::setw(9);
            out<<"avg,";
            out<<std::right<<std::setw(3);
            out<<key;
            out<<std::left<<std::setw(9);
            out<<"std,";
        }
        out<<std::endl;

        out<<std::left<<std::setw(11);
        out<<repetitions<<",";
        for (int key : iter_pastkeys)
        {
            out<<std::left<<std::setw(11);
            std::stringstream stream0;
            stream0 << std::fixed << std::setprecision(7) << iter_avgs[key];
            out<<stream0.str()<<",";
            out<<std::left<<std::setw(11);
            std::stringstream stream1;
            stream1 << std::fixed << std::setprecision(7) << iter_stds[key];
            out<<stream1.str()<<",";
        }
        out<<std::endl;
        out.close();   
    }

    // 2.2 If file existed, rewrite file with old and new info
    else
    {
        int N_prev_columns = 0;
        std::string line_0 = previous_lines[0];
        for (int j = 0; j<line_0.size(); j++){
            if (line_0[j] == ',') N_prev_columns += 1;}
        int current_keys = iter_pastkeys.size();
        int n_prev_keys = (N_prev_columns-1)/2;
        int delta = current_keys - n_prev_keys;
        
        // If file had less columns, write new one with old data + new columns
        if (delta > 0)
        {
            std::ofstream out(filename);
            for (int i = 0; i < previous_lines.size(); i++)
            {
                std::string line_i = previous_lines[i];
                // complete labels
                if (i==0){
                    out<<line_i;
                    for (int x = n_prev_keys; x<current_keys; x++)
                    {
                        int key = iter_pastkeys[x];
                        out<<std::right<<std::setw(3);
                        out<<key;
                        out<<std::left<<std::setw(9);
                        out<<"avg,";
                        out<<std::right<<std::setw(3);
                        out<<key;
                        out<<std::left<<std::setw(9);
                        out<<"std,";
                    }
                    out<<std::endl;
                }
                // add columns, i.e. zeros at every line
                else
                {
                    out<<line_i;
                    std::vector<double> zeroes (delta*2, 0.);
                    for (double zero : zeroes)
                    {
                        out<<std::left<<std::setw(11);
                        out<<zero<<",";
                    }
                    out<<std::endl;
                }       
            }
            out.close();
        }

        // append new lines
        std::ofstream out(filename, std::ios::app);
        out<<std::left<<std::setw(11);
        out<<repetitions<<",";
        for (int i = 0; i<iter_pastkeys.size(); i++)
        {
            int key = iter_pastkeys[i];
            out<<std::left<<std::setw(11);
            std::stringstream stream0;
            stream0 << std::fixed << std::setprecision(7) << iter_avgs[key];
            out<<stream0.str()<<",";
            out<<std::left<<std::setw(11);
            std::stringstream stream1;
            stream1 << std::fixed << std::setprecision(7) << iter_stds[key];
            out<<stream1.str()<<",";
        }
        if (delta<0)
        {
            std::vector<double> zeroes (-delta*2, 0.);
            print_vector(zeroes);
            for (double zero : zeroes){
                out<<std::left<<std::setw(11);
                out<<zero<<",";
            }
        }
        out<<std::endl;
        out.close();
    }
}


auto finish = high_resolution_clock::now();
double duration = duration_cast<microseconds>(finish - beginning).count();
std::cout<<"\n=============================================================\n";
std::cout << "\nProgram took in total: "
        << duration/pow(10,6) << " seconds\n" << std::endl;

std::cout << "Parameters used:\n";
std::cout << "Dim = "<< dim << ", Rho = "<< Rho << std::endl;
}








///////////////////////////////////////////////////////////////////////////////
// USEFUL FUNCTIONS' DEFINITIONS
//////////////////////////////////////////////////////////////////////////////
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
                  int N)
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