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
        // Set filename
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

        std::string filename_base = std::string(homeDir) 
                    + "/MSci_Schwarzschild_Causets/data/HRVs_connectivity/"
                    + "M=" + mass_str
                    + "_Rho=" + rho_str 
                    + "_Card=" + std::to_string(cards[0])
                    + "_r=" + radius_str
                    + "_hollow=" + hollow_str
                    + "_dur=" + dur_str;
                                
        // Get ith number of file (how many existed before)
        int i = 0;
        bool hasopened = true;
        while (hasopened)
        {
            std::string filename_i =filename_base + std::to_string(i) + ".txt";

            std::vector<std::string> previous_lines;
            std::string line;
            std::ifstream prev_file(filename_i);
            while (std::getline(prev_file, line)) 
            {
                previous_lines.push_back(line);
            }
            prev_file.close();

            // If file didn't exist, break -> you have found how many existed
            if (previous_lines.size()==0){
                hasopened = false;
                break;
            }
            else{
                i+=1;
            }
        }
        std::string filename_new = filename_base + std::to_string(i) + ".txt";

        
    
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

        // Count HRVs and save in file
        auto countstart = high_resolution_clock::now();
        double t_f = 0; 
        C.save_molecules_only(filename_new.c_str(), 0, 2*mass, "HRVs");
        auto countend = high_resolution_clock::now();
        double durationlinks = duration_cast<microseconds>(
                                countend - countstart).count();
        std::cout << "Time taken in count_HRVs for N = "
        << C._size << ": " << durationlinks/pow(10,6) << " seconds"
        << std::endl;  

        std::cout << "Saving Iteration "<< iteration <<
                    "to the file: " << filename_new << std::endl;           
    } 
}

auto finish = high_resolution_clock::now();
double duration = duration_cast<microseconds>(finish - beginning).count();
std::cout<<"\n=============================================================\n";
std::cout<<"===============================================================\n";
std::cout << "\nProgram took in total: "
        << duration/pow(10,6) << " seconds\n" << std::endl;

std::cout << "Parameters used:\n";
}
