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

double x = xi_0(4,2);
double y = chi_k(4,2);
std::cout << "xi=" << x << std::endl;
std::cout << "chi=" << y << std::endl;


// Create causet in Schwarzschild spacetime
///////////////////////////////////////////////////////////////////////////////
                
std::vector<int> dims = {4}; 
std::vector<int> cards = {1000};
double mass = 0.25;
int N_reps = 1;


// Sprinkling Parameters
///////////////////////////////////////////
bool poisson = true;
bool make_matrix = true;
bool special = false;
bool use_transitivity = false;
bool make_sets = false;
bool make_links = false; 
const char* sets_type = "all"; 
const char* name = "cylinder";

// Shape parameters
double radius = 1;
double height = 2;
///////////////////////////////////////////

// Begin program
auto beginning = high_resolution_clock::now();

std::cout<<"\n\n============= Sprinkling into "<<name<<" ==================\n";

int iteration = 0;
for (auto dim: dims)
{
    auto start = high_resolution_clock::now();

    std::cout << "Dim = " << dim << std::endl;
    
    for (auto card : cards)
    {
        for (int rep=0; rep<N_reps; rep++)
        {
                auto repstart = high_resolution_clock::now();
                // Set up shape
                std::vector<double> center = {0.0,0.0,0.0,0.0};
                CoordinateShape shape(dim,name,center,radius,height);
                // Set up spacetime
                Spacetime S = Spacetime();
                S.BlackHoleSpacetime(dim,mass);
                // Sprinkle the causet
                SprinkledCauset C(card, S, shape, poisson,
                                make_matrix, special, use_transitivity,
                                make_sets, make_links,sets_type);

                
                std::cout << "Dim="<< dim <<", "<<(rep+1)<<"/"<<N_reps<<"\n";

                //Timing rep
                auto repend = high_resolution_clock::now();
                double duration = duration_cast<microseconds>(repend - repstart).count();
                std::cout << "Time taken generating for N = " << C._size
                << ": " << duration/pow(10,6) << " seconds" << std::endl;
        }
    
        auto mid = high_resolution_clock::now();
        double duration = duration_cast<microseconds>(mid - start).count();
        
        std::cout << "Average time taken for generating "<< N_reps
                << " causets with N = " << card << ":\n"  
                << duration/pow(10,6)/N_reps
                << " seconds\n" << std::endl;
    }   
}


auto finish = high_resolution_clock::now();
double duration = duration_cast<microseconds>(finish - beginning).count();
std::cout<<"===============================================================\n";
std::cout << "\nProgram took in total: "
        << duration/pow(10,6) << " seconds\n" << std::endl;



// end of main   
}