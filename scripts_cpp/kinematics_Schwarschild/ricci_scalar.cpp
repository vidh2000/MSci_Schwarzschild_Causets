#include <algorithm>
#include <cmath>
#include <chrono>
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
#include <unordered_set>

#include "../causets_cpp/sprinkledcauset.h"
#include "../causets_cpp/shapes.h"
#include "../causets_cpp/spacetimes.h"

#include "../causets_cpp/functions.h"
#include "../causets_cpp/vecfunctions.h"

#include "../causets_cpp/kinematics_functions.h"

//#include <boost/range/combine.hpp>
//#include <boost/math/special_functions/gamma.hpp>
#include <omp.h>

// $HOME var get
#include <unistd.h>
#include <sys/types.h>
//#include <pwd.h>

using namespace std::chrono;
using namespace boost::math;





// Defining Ricci scalar (components) and proper time extension of the interval
///////////////////////////////////////////////////////////////////////////////





int main(){


// Create causet in Schwarzschild spacetime
///////////////////////////////////////////////////////////////////////////////
                
std::vector<int> dims = {4}; 
std::vector<int> cards = {5000};
int min_size = 25;  //Minimal size of the interval (min # of elements in it)
int max_size = 0;   //if 0 -> ==_size of the causet
double mass = 0.25;
int N_reps = 50;
int N_intervals = 10; // Number of intervals per causet realisation


// Sprinkling Parameters
///////////////////////////////////////////
bool poisson = true;
bool make_matrix = true;
bool special = false;
bool use_transitivity = false;
bool make_sets = false;
bool make_links = false; 
const char* sets_type = "future";
const char* name = "cylinder";

// Shape parameters
double radius = 4*mass;
double height = 1; 
double hollow = 0.2;
///////////////////////////////////////////

// Begin program
auto beginning = high_resolution_clock::now();

std::cout<<"\n\n============= Sprinkling into "<<name<<" ==================\n";
int iteration = 0;
for (auto dim: dims)
{
    auto start = high_resolution_clock::now();

    for (auto card : cards)
    {
        /*Array of pairs for all intervals (from different causets):
              item: <N_chains vector of the interval, r_avg of the interval> */
        std::vector<std::pair<std::vector<double>,double>> C_k_arr;


        // Generate causets and find intervals and number of k-chains...
        ////////////////////////////////////////////////////////////////////
        int rep=0;
        while (rep<N_reps)
        {
            std::cout << "Dim="<< dim <<", "<<(rep+1)<<"/"<<N_reps<<"\n";
            auto repstart = high_resolution_clock::now();
            // Set up shape
            std::vector<double> center = {0.0,0.0,0.0,0.0};
            CoordinateShape shape(dim,name,center,radius,height,hollow);
            // Set up spacetime
            Spacetime S = Spacetime();
            S.BlackHoleSpacetime(dim,mass);
            // Sprinkle the causet
            SprinkledCauset C(card, S, shape, poisson,
                            make_matrix, special, use_transitivity,
                            make_sets, make_links,sets_type);

            double nRels = sumMatrix(C._CMatrix);
            std::cout<<"Number of relations:" << nRels << std::endl;
            
            // Get array of "N_chains" for chain-sizes 1...4.
            std::vector<std::pair<std::vector<double>,double>> nchains_arr;
            try {
                nchains_arr = C.get_Nchains_inInterval(N_intervals,
                                            min_size, 4, max_size);
            } catch (std::runtime_error& e) {
                continue;
            }
            
            // Store Nchain_k vectors and r_avg values 
            for (auto item : nchains_arr) {
                print(item.first);
                C_k_arr.push_back(item);
            }

            //Timing rep
            auto repend = high_resolution_clock::now();
            double duration = duration_cast<microseconds>(repend - repstart).count();
            std::cout << "Rep finished, N = " << C._size
            << ": " << duration/pow(10,6) << " seconds" << std::endl;

            // Add to counter of repetitions
            rep++;

        }
    
        auto mid = high_resolution_clock::now();
        double duration = duration_cast<microseconds>(mid - start).count();
        
        std::cout << "-------------------------------------\n";
        std::cout << "<t> for "<< N_reps
                << " causets in D=" << dim << " with N = " << card << ": "  
                << duration/pow(10,6)/N_reps
                << " seconds\n" << std::endl;



        // Find the R_scalar and R_tensor and T (proper time.. if)
        //////////////////////////////////////////////////////////////////
        print("Finding Ricci scalar and tensor 00th component"); //also find T?

        // Calculate the density - assume we always sprinkle into cylinder
        double density;
        double volume;
        if (dim==3) {
            volume = pi*radius*radius*height;
        }
        else if (dim==4) {
            volume = 4/3*pi*radius*radius*radius*height;
        }
        else {
            print("ERROR: Dimension must be 3 or 4!");
            throw std::runtime_error("");
        }
        density = card/volume;

        // Find values of causet kinematics-related variables
        std::vector<double> R_scalar_arr;
        std::vector<double> R_tensor00_arr;

        for (auto item : C_k_arr) {
            // item == <array of Nchains_k, r_avg>
            R_scalar_arr.push_back(R_RSS(dim, item.first, density));
            R_tensor00_arr.push_back(R_00(dim, item.first, density)); 
        } 
        
        print("\nxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n");
        std::cout << "Density = " << density << std::endl;
        print("Using units of discreteness scale");

        double R_RSS_scalar_avg = mymean(R_scalar_arr)/std::sqrt(density);
        double R_RSS_scalar_std = mystd(R_scalar_arr)/std::sqrt(density);
        
        double R_RSS_tensor00_avg = mymean(R_tensor00_arr)/std::sqrt(density);
        double R_RSS_tensor00_std = mystd(R_tensor00_arr)/std::sqrt(density);
        

        std::cout << "Dim="<<dim<<". Card="<<card << std::endl;
        

        std::cout << "Rscalar = " << R_RSS_scalar_avg
            << "+-" << R_RSS_scalar_std << std::endl;
        
        std::cout << "R_tensor (00th component)  = "
            << R_RSS_tensor00_avg << "+-"
            << R_RSS_tensor00_std << std::endl;

        // Finished for that Dimension and Cardinality value
    }   
}


auto finish = high_resolution_clock::now();
double duration = duration_cast<microseconds>(finish - beginning).count();
std::cout<<"===============================================================\n";
std::cout << "\nProgram took in total: "
        << duration/pow(10,6) << " seconds\n" << std::endl;




// end of main   
}