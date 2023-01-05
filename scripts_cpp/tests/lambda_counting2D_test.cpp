
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

using namespace std;
using namespace std::chrono;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


////////////////////////////////
int dim = 2;
double t_f = 0;
double r_S = 2;
std::vector<int> cards = {25, 50, 75, 100, 200, 300};

// Cube Shape Parameters
double radius = 5;
double myduration = 5;
vector<double> loop_edges = {2.5}; //edges of cubes to loop over

// Sprinkle Parameters
bool poisson          = false;
bool make_matrix      = false;
bool special          = false;
bool use_transitivity = true;
bool make_sets        = true;
bool make_links       = false;
const char* sets_type = "all with links";



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

int main()
{
    auto beginning = high_resolution_clock::now();

    for (int card : cards)
    {
        std::cout<<"\nCARD "<<card<<" \n";
        for (double edge : loop_edges)
        {
            std::cout<<"EDGE "<<edge<<" \n";

            // Set up shape
            vector<double> center (dim, edge/2); 
            center[0] -= edge;
            //center[0] -= edge/2;
            std::vector<const char*> shapes = {"cylinder", "cube"};
            int shapeindex = 1; (dim==2)? 0 : 1;
            CoordinateShape Shape(dim,shapes[shapeindex],center,
                                radius, myduration, edge);

            // Set up spacetime
            Spacetime S = Spacetime();
            S.BlackHoleSpacetime(dim, r_S/2);

            // Sprinkle the causet
            SprinkledCauset C(card, S, Shape, poisson,
                                make_matrix, special, use_transitivity,
                                make_sets, make_links,sets_type);

            //Save Lambdas
            std::stringstream estream;
            estream << std::fixed << std::setprecision(2) << edge;
            std::stringstream rstream;
            rstream << std::fixed << std::setprecision(2) << radius;
            std::string redge_s = (dim > 0)? estream.str() : rstream.str();

            std::string path_file_str = "../../data/blackhole_and_lambdas"
                                            + std::to_string(dim)
                                            + "D_N" + std::to_string(card)
                                            + "_redge" + redge_s + ".txt";
            const char* path_file = path_file_str.c_str();

            C.save_lambdas(path_file, "sets", t_f, r_S);
        }
    }

auto finish = high_resolution_clock::now();
double duration = duration_cast<microseconds>(finish - beginning).count();
std::cout << "\nProgram took in total: "
        << duration/pow(10,6) << " seconds\n" << std::endl;

}