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


///////////////////////////////////////////////////////////////////////////////
// USEFUL FUNCTIONS (defined below everything)
///////////////////////////////////////////////////////////////////////////////

template <typename num>
std::vector<std::map<int, double>> avg_distr(
                                std::map<int, std::vector<num>> &all_results);

template <typename num>
void update_distr(std::map<int, std::vector<num>> &all_results, 
                  std::map<int, num> &newresults);



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//        PARAMETERS are inputted via COMMAND LINE in BASH script
//              - param 1 = mass (double)
//              - param 2 = Rho(double)
//              - param 3 = N_reps (int)
//___________________________________________________________________________//
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
double R = 2*mass+4*scale;
std::cout << "R done"<< std::endl;
double r = 2*mass-4*scale;
std::cout << "r done"<< std::endl;
double T = 6*scale;
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
auto beginning = high_resolution_clock::now();

std::cout<<"\n\n============ Sprinkling into "<<name<<" ===================\n";
std::cout << "Doing CMatrix and inferring future links from it\n \n";

// Variables for storage of information from each iteration
std::map<int, std::vector<double>> all_HRVs_results;

int iteration = 0;
for (auto && tup : boost::combine(cards, radii, hollow_vals,
                                    masses, durations, repetitions_arr))
{
    iteration++;
    auto start = high_resolution_clock::now();

    // Store N_links for each repetition in a single iteration
    std::map<int, std::vector<double>> all_iter_HRVs_results;
    std::vector<int> HRV_keys = {-3, -2, -1, 0, 1};
    
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

            // Count HRVs and update results
            auto countstart = high_resolution_clock::now();
            double t_f = 0;
            std::map<int, double> HRVs_distr = C.count_HRVs(t_f,2*mass);
            update_distr(all_iter_HRVs_results, HRVs_distr);

            //Timing link counting
            auto countend = high_resolution_clock::now();
            double durationlinks = duration_cast<microseconds>(
                                    countend - countstart).count();
            std::cout << "Time taken in count_HRVs for N = "
            << C._size << ": " << durationlinks/pow(10,6) << " seconds"
            << std::endl;     
    }

    std::vector<std::map<int, double>> iter_results = avg_distr(
                                                    all_iter_HRVs_results);
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
                            + "/MSci_Schwarzschild_Causets/data/HRVs/"
                            + "M=" + mass_str
                            + "_Rho=" + rho_str 
                            + "_Card=" + std::to_string(cards[0])
                            + "_r=" + radius_str
                            + "_hollow=" + hollow_str
                            + "_dur=" + dur_str
                            + ".txt";
    std::cout<<"\n==========================================================\n";
    std::cout << "Saving Iteration "<< iteration <<
                "to the file: " << filename << std::endl;  
    

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

    // 2 If file didn't exist, write for first time
    if (previous_lines.size()==0)
    {
        std::ofstream out(filename);
        // First line: labels
        out<<std::left<<std::setw(11);
        out<<"Nreps,";
        out<<std::left<<std::setw(11);
        out<<"outer_avg,";
        out<<std::left<<std::setw(11);
        out<<"outer_std,";
        out<<std::left<<std::setw(11);
        out<<"inner_avg,";
        out<<std::left<<std::setw(11);
        out<<"inner_std,";
        out<<std::left<<std::setw(11);
        out<<"mintm_avg,";
        out<<std::left<<std::setw(11);
        out<<"mintm_std,";
        out<<std::left<<std::setw(11);
        out<<"open_avg,";
        out<<std::left<<std::setw(11);
        out<<"open_std,";
        out<<std::left<<std::setw(11);
        out<<"close_avg,";
        out<<std::left<<std::setw(11);
        out<<"close_std,";
        out<<std::endl;
        out.close();   
    }

    // 3 Append new info
    std::ofstream out(filename, std::ios::app);
    out<<std::left<<std::setw(11);
    out<<repetitions<<",";
    for (int key : HRV_keys)
    {
        out<<std::left<<std::setw(10);
        std::stringstream stream0;
        stream0 << std::fixed << std::setprecision(6) << iter_avgs[key];
        out<<stream0.str()<<",";
        out<<std::left<<std::setw(10);
        std::stringstream stream1;
        stream1 << std::fixed << std::setprecision(6) << iter_stds[key];
        out<<stream1.str()<<",";
    }
    out<<std::endl;
    out.close();
}


auto finish = high_resolution_clock::now();
double duration = duration_cast<microseconds>(finish - beginning).count();
std::cout<<"\n=============================================================\n";
std::cout<<"===============================================================\n";
std::cout << "\nProgram took in total: "
        << duration/pow(10,6) << " seconds\n" << std::endl;

std::cout << "Parameters used:\n";
}












//////////////////////////////////////////////////////////////////////////////
// USEFUL FUNCTIONS' DEFINITIONS
//////////////////////////////////////////////////////////////////////////////



/**
 * @brief Add new results to previous, whereeach results is given by a map
 * <int, vector<int>>, such as for HRVs.
 * 
 * @param all_results 
 * @param newresults 
 */
template <typename num>
void update_distr(std::map<int, std::vector<num>> &all_results, 
                  std::map<int, num> &newresults)
{                
    //Update past results with new ones
    for (int key : {-3, -2, -1, 0, 1})
    {
        all_results[key].push_back(newresults[key]);
    }
}


/**
 * @brief get {avgs, stds}, where avgs and stds are maps int->double, where
 * int is the key (size of HRVs) and double is the avg/std of the results (the 
 * counts of the HRVs of that size).
 * 
 * @param all_results 
 */
template <typename num>
std::vector<std::map<int, double>> avg_distr(
                                std::map<int, std::vector<num>> &all_results)
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