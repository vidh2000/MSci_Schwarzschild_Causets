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

#include "../scripts_cpp/causets_cpp/sprinkledcauset.h"
#include "../scripts_cpp/causets_cpp/shapes.h"
#include "../scripts_cpp/causets_cpp/spacetimes.h"

#include "../scripts_cpp/causets_cpp/functions.h"
#include "../scripts_cpp/causets_cpp/vecfunctions.h"

using namespace std::chrono;

int card = 200;
int dim = 2;
std::vector<double> center (dim, 10);
double radius = 2.0;
double myduration = 10;
double edge = 20;
std::vector<double> edges = {1,2,3,4};

// Parameters
bool poisson = false;
bool make_matrix = false;
bool special = false;
bool use_transitivity = true;
bool make_sets = true;
bool make_links = false;
const char* sets_type = "all with links";

int main(){

auto start = high_resolution_clock::now();

std::vector<const char*> names = {"cube"};

for (const char* name : names)
{
    std::cout<<"\n\n============= USING "<<name<<" ====================\n";

    CoordinateShape shape(dim,name,center,radius,edge,edges,0.0,myduration);
    Spacetime S = Spacetime();
    //S.FlatSpacetime(dim);
    S.BlackHoleSpacetime(dim=2);
    SprinkledCauset C(card, S, shape, poisson,
                        make_matrix, special, use_transitivity,
                        make_sets, make_links,sets_type);
    std::cout << "Generate the causet... Saving ->" << std::endl;
    const char* path_file = "../data/blackhole2D_N200_edge20.txt";
    C.save_causet(path_file);
}

auto stop = high_resolution_clock::now();
    double duration = duration_cast<microseconds>(stop - start).count();
    std::cout << "Time taken: "
            << duration/pow(10,6) << " seconds" << std::endl;

}