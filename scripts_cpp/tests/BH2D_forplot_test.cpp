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

#include "../causets_cpp/embeddedcauset.h"
#include "../causets_cpp/sprinkledcauset.h"
#include "../causets_cpp/shapes.h"
#include "../causets_cpp/spacetimes.h"
#include "../causets_cpp/functions.h"
#include "../causets_cpp/vecfunctions.h"

using std::cout;
using std::endl;
using std::vector;
using namespace std::chrono;

////////////////////////////////
std::vector<int> cards = {25, 50, 75, 100, 150, 250, 500};//, 150, 200, 500};
int dim = 2;
double mass = 1;
vector<std::string> metrics = {"EF(original)", "EF(uv)", "S"} ;

// Shape Parameters
const char* name = "cube";
double radius = 3.0;
double myduration = 3;
double edge = 2;
std::vector<double> center (dim, edge/2);
std::vector<double> edges = {1,2,3,4};

// Sprinkle Parameters
bool poisson = false;
bool make_matrix = false;
bool special = false;
bool use_transitivity = true;
bool make_sets = true;
bool make_links = false;
const char* sets_type = "all with links";

int main()
{
    // Set string version of edge and radius to 2 precision for saving names
    std::stringstream estream;
    estream << std::fixed << std::setprecision(2) << edge;
    std::string edge_s = estream.str();
    std::stringstream rstream;
    rstream << std::fixed << std::setprecision(2) << radius;
    std::string radius_s = rstream.str();

    //Start
    auto start = high_resolution_clock::now();

    std::cout<<"\n\n============= USING "<<name<<" ====================\n";
    for (std::string metric : metrics)
    {
        for (int card : cards)
        {
            cout<<"Cardinality = "<<card<<endl;
            CoordinateShape shape (dim,name,center,radius,myduration,
                                   edge,edges,0.0);
            Spacetime S;
            S.BlackHoleSpacetime(dim, mass, metric);
            SprinkledCauset C(card, S, shape, poisson,
                                make_matrix, special, use_transitivity,
                                make_sets, make_links,sets_type);
            std::cout << "Generated the causet... Saving ->" << std::endl;

            if (metric == "EF(uv)"){
            std::string path_file_str = "../../data/blackhole_EFv_"
                                        + std::to_string(dim)
                                        + "D_N" + std::to_string(card)
                                        + "_redge_";
            path_file_str += (std::strcmp(name, "cube")==0)? edge_s : radius_s;
            path_file_str += ".txt";
            const char* path_file = path_file_str.c_str();
            C.save_causet(path_file);}

            else if (metric == "S"){
            std::string path_file_str = "../../data/blackhole_S_"
                                        + std::to_string(dim)
                                        + "D_N"+std::to_string(card)
                                        + "_redge_";
            path_file_str += (std::strcmp(name, "cube")==0)? edge_s : radius_s;
            path_file_str += ".txt";
            const char* path_file = path_file_str.c_str();
            C.save_causet(path_file);}

            else{
            std::string path_file_str = "../../data/blackhole"
                                        + std::to_string(dim)
                                        + "D_N" + std::to_string(card)
                                        + "_redge_";
            path_file_str += (std::strcmp(name, "cube")==0)? edge_s : radius_s;
            path_file_str += ".txt";
            const char* path_file = path_file_str.c_str();
            C.save_causet(path_file);}
        }
    }

    auto stop = high_resolution_clock::now();
    double duration = duration_cast<microseconds>(stop - start).count();
    std::cout << "Time taken: "<< duration/pow(10,6) << " s" << std::endl;
}