
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
vector<int> dims = {3, 4};
double t_f = 0;
double r_S = 2;
std::vector<int> cards = {1000, 2000, 3000};

// Cube Shape Parameters
double radius = 4;
double myduration = 4;
vector<double> loop_edges = {4}; //edges of cubes to loop over 
bool centre_cube_in_horizon = false;//(spatially center 2D cube on horizon)

// Sprinkle Parameters
bool poisson          = true;
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

    for (int dim : dims)
    { 
    for (int card : cards)
    {
        for (double edge : loop_edges)
        {
            // Set up shape
            std::vector<const char*> shapes = {"cylinder", "cube"};
            vector<double> centre (dim, 0);
            int shapeindex = 0;
            if (dim == 2)
            {
                shapeindex += 1; //pick cube
                centre[0] -= edge/2;
                if (centre_cube_in_horizon)
                    centre[1] = r_S;
                std::cout<<" \nDIM "<<dim
                         <<" ;CARD "<<card
                         <<" ;EDGE "<<edge
                         <<"\n";
            }
            else
            {
                centre[0] -= myduration/2;
                std::cout<<" \nDIM "<<dim
                         <<" ;CARD "<<card
                         <<" ;RADIUS  "<<radius
                         <<" ;DURATION "<<myduration
                         <<"\n";
            }
            CoordinateShape Shape(dim,shapes[shapeindex],centre,
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
            std::string redge_s = (shapeindex==1)?estream.str():rstream.str();

            std::string path_file_str = "../../../data/data_for_plotting/blackhole_and_HRVs"
                                            + std::to_string(dim)
                                            + "D_N" + std::to_string(card)
                                            + "_redge" + redge_s;
            if (centre_cube_in_horizon)
                path_file_str += "_horizon_centred";
            path_file_str +=  ".txt";
            const char* path_file = path_file_str.c_str();

            std::cout<<"Counting links"<<std::endl;
            C.save_molecules(path_file, "sets", t_f, r_S, "HRVs");
        }
    }
    }

auto finish = high_resolution_clock::now();
double duration = duration_cast<microseconds>(finish - beginning).count();
std::cout << "\nProgram took in total: "
        << duration/pow(10,6) << " seconds\n" << std::endl;

}