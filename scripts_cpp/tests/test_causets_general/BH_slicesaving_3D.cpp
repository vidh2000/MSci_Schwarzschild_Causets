
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

#include "../../causets_cpp/sprinkledcauset.h"
#include "../../causets_cpp/shapes.h"
#include "../../causets_cpp/spacetimes.h"

#include "../../causets_cpp/functions.h"
#include "../../causets_cpp/vecfunctions.h"

#include <boost/range/combine.hpp>
#include <omp.h>

// $HOME var get
#include <unistd.h>
#include <sys/types.h>

using namespace std;
using namespace std::chrono;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


////////////////////////////////
vector<int> dims = {3};
double t_f = 0;
double r_S = 2;
double deltaphi = 2*3.1415 / (2*1440);
std::vector<int> cards = {50, 60, 75, 90, 100, 125}; 

// Sprinkle Parameters
bool poisson          = true;
bool make_matrix      = true;
bool special          = false;
bool use_transitivity = false;
bool make_sets        = false;
bool make_links       = false;
const char* sets_type = "all with links";



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

int main()
{
    auto beginning = high_resolution_clock::now();

    for (int dim : dims)
    { 
    for (int card : cards)
    {
        double mass = r_S /2.;
        double R = 2*r_S;
        double T = R*1.;
        double h = 0;
        std::cout << "\nDim = " << dim << std::endl;
        std::cout << "Card = " << card << std::endl;
        
        // Set up shape
        const char* shape = "cylinder";
        vector<double> centre (dim, 0);
        centre[0] -= T/2;
        CoordinateShape Shape(dim,shape,centre,R, T, h);
        Shape.limitphi(deltaphi);

        // Set up spacetime
        Spacetime S = Spacetime();
        S.BlackHoleSpacetime(dim, r_S/2);

        // Sprinkle the causet
        SprinkledCauset C(card, S, Shape, poisson,
                            make_matrix, special, use_transitivity,
                            make_sets, make_links,sets_type);

        //Save Lambdas
        const char* homeDir = getenv("HOME");
        std::stringstream rstream;
        rstream << std::fixed << std::setprecision(2) << R;
        std::string redge_s = rstream.str();
        // std::stringstream hstream;
        // hstream << std::fixed << std::setprecision(2) << h;
        // std::string h_s = hstream.str();

        std::string filename = "../../../data/data_for_plotting/blackhole_slice"
                            + std::to_string(dim)
                            + "D_N" + std::to_string(card)
                            + "_redge" + redge_s;
        filename +=  ".txt";
        const char* path_file = filename.c_str();
        std::cout<<"Ready to save in\n" << path_file << std::endl;
        C.save_causet(path_file, "sets");
    }
    }

    auto finish = high_resolution_clock::now();
    double duration = duration_cast<microseconds>(finish - beginning).count();
    std::cout << "\nProgram took in total: "
            << duration/pow(10,6) << " seconds\n" << std::endl;
}
