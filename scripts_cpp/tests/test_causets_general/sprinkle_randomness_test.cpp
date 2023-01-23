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

using namespace std::chrono;//can't use duration

//////////////////////////////////////////////////////////////////////
// for a complete test of randomness in coordinates, it is advisable
// to comment out this->sort_coords and this->make_attrs in 
// SprinkledCauset (lines 86-88)
// THEN PUT THEM ON AGAIN!!!!!
/////////////////////////////////////////////////////////////////////

// Sprinkled causet parameters
int card = 50;
int repetitions = 50;
int dim = 3;
std::vector<double> center (dim, 0.0);
double radius = 4.0;
double myduration = 10;
double edge = 1.5;
double hollow = 0;
std::vector<double> edges = {1,2,3,4};


// Sprinkling Parameters
bool poisson = true;
bool make_matrix = true;
bool special = false;
bool use_transitivity = false;
bool make_sets = false;
bool make_links = false; 
const char* sets_type = "future"; 
const char* name = "cylinder";

//Analysis Parameters (printing it?)
bool want_coords = false;


int main(){
    
    auto start = high_resolution_clock::now();
    //std::cout << "Starting building shape..." << std::endl;
    std::vector<const char*> names = {"cylinder"};// "ball", "cylinder", "cube", "cuboid"};
    edges.resize(dim);
    //card = (want_coords || want_matrix)? 10 : card;
    for (const char* name : names)
    {

        CoordinateShape shape(dim,name,center,radius,myduration,
                            hollow,edge,edges);

        std::cout<<"\n\n============= USING "<<name<<" ====================\n";
        // std::cout << "What are the just-Shape's parameters?\n";
        // for (auto const& p : shape._params)
        //     {std::cout << "-- " << p.first << "=" << p.second << '\n';}
        

        Spacetime S = Spacetime();
        S.FlatSpacetime(dim);

        std::vector<int> cards = {};
        std::vector <std::vector< std::vector<double> > > all_coords = {};
        for (int r = 0; r<repetitions; r++)
        {
            SprinkledCauset C(card, S, shape, poisson,
                                make_matrix, special, use_transitivity,
                                make_sets, make_links,sets_type);
            cards.push_back(C._size);

            if (want_coords){
            std::cout << "\nCoordinates:\n";
            print_vector(C._coords);}
            all_coords.push_back(C._coords);
        }

        // ANALYSE POISSON DISTRIBUTION OF CARDINALITIES
        double avg_card = mymean(cards);
        double accum = 0.0;
        std::for_each(std::begin(cards), std::end(cards),
                        [&](double N){accum += (N - avg_card) * (N - avg_card);}
                        );
        double std_card = sqrt(accum / cards.size());
        std::cout<<"The Cardinalities are :";
        print_vector(cards);
        std::cout<<"With average : "<<avg_card<<std::endl;
        std::cout<<"With std     : "<<std_card<<std::endl;

        //CHECK COORDS VARY
        std::vector<std::vector<double>> avg_coords(vecmin(cards), std::vector<double>(dim,0));
        std::vector<std::vector<double>> std_coords(vecmin(cards), std::vector<double>(dim,0));
        for (int i = 0; i < vecmin(cards); i++) 
        {
            // first average
            for (int mu = 0; mu < dim; mu++)
            {
                for (int r = 0; r<repetitions; r++)
                {
                    avg_coords[i][mu] += all_coords[r][i][mu];
                }
                avg_coords[i][mu] /= repetitions;
            }

            // then std
            for (int mu = 0; mu < dim; mu++)
            {
                for (int r = 0; r<repetitions; r++)
                {
                    std_coords[i][mu] += (all_coords[r][i][mu]-avg_coords[i][mu])*
                                         (all_coords[r][i][mu]-avg_coords[i][mu]);
                }
                std_coords[i][mu] = std::sqrt(std_coords[i][mu]/repetitions);
            }
        }
        std::cout<<"\nThe Coordinates Averages are :\n";
        print_vector(avg_coords);
        std::cout<<"\nThe Coordinates Stds are :\n";
        print_vector(std_coords);
    }
    auto stop = high_resolution_clock::now();
    double duration = duration_cast<microseconds>(stop - start).count();
    std::cout << "Time taken for "<<repetitions<<"repetitions, N= " << card 
            << ", dim = "<< dim << ":\n"
            << duration/pow(10,6) << " seconds" << std::endl;

    std::cout << "\n =========This file works!==========\n" << std::endl;


    return 0;
};