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

#include "causets_cpp/sprinkledcauset.h"
#include "causets_cpp/shapes.h"
#include "causets_cpp/spacetimes.h"

#include "causets_cpp/functions.h"
#include "causets_cpp/vecfunctions.h"

using namespace std::chrono;
/**
 * @brief   To run this file, you must:
 *          - Go to to the folder in which it is located
 *          - Type into the command line:
 *              cd scripts_cpp
 *              g++ -g causets_cpp/causet.cpp causets_cpp/shapes.cpp causets_cpp/spacetimes.cpp causets_cpp/embeddedcauset.cpp causets_cpp/sprinkledcauset.cpp sprinkle.cpp -std=c++17 -o sprinkle -O2
 *          - This creates the executable sprinkle.exe
 *          - Run sprinkle.exe by typing into command line:
 *              ./sprinkle
 * 
 *          NOTE: running in Windows cmd does not print no matetr what
 * 
 */
// gives 0 but doesn't create a new entry std::cout << "sprinkle.cpp 'shape' radius= " << shape._params.find("radius")->second << std::endl;
// gives 0 but creates a new entry std::cout << "sprinkle.cpp 'shape' radius= " << shape._params["radius"]<< std::endl;
   

// Sprinkled causet parameters
int card = 5000;
int dim = 4;
std::vector<double> center (dim, 0.0);
double radius = 4.0;
double myduration = 10;
double edge = 1.5;
std::vector<double> edges = {1,2,3,4};

//Analysis Parameters
bool want_matrix = false;
bool want_coords = false;

bool poisson = false;
bool make_matrix = true;
bool special = true;
bool use_transitivity = false;
bool make_sets = true;
bool make_links = false;

int main(){
    auto start = high_resolution_clock::now();
    //std::cout << "Starting building shape..." << std::endl;
    std::vector<const char*> names = {"bicone"};// "ball", "cylinder", "cube", "cuboid"};
    edges.resize(dim);
    card = (want_coords || want_matrix)? 10 : card;
    for (const char* name : names)
    {

        std::cout<<"\n\n============= USING "<<name<<" ====================\n";
        std::cout << "What are the just-Shape's parameters?\n";
        CoordinateShape shape(dim,name,center,radius,edge,edges,0.0,myduration);
        for (auto const& p : shape._params)
            {std::cout << "-- " << p.first << "=" << p.second << '\n';}
        

        Spacetime S = Spacetime();
        S.FlatSpacetime(dim);
        SprinkledCauset C(card, S, shape, poisson,
                          make_matrix, special, use_transitivity,
                          make_sets, make_links);
        C.make_futures();
        
        //PARAMETERS
        //std::cout << "\nWhat are the Causet's Shape's parameters at the end?\n";
        //for (auto const& p : C._shape._params){
        //    std::cout << "-- " << p.first << '=' << p.second << '\n';}
        if ((strcmp(shape._name, "ball")==0)      ||
             (strcmp(shape._name, "cylinder")==0) ||
             (strcmp(shape._name, "diamond")==0)  ||
             (strcmp(shape._name, "bicone")==0))
        {
            std::cout<<"Radius obtained from Parameter: "
                 <<C._shape.Parameter("radius")<<std::endl;
        }
        // std::cout<<"-- cardinality=" << C._size << "\n";
        // std::cout<<"-- dim=" << C.spacetime_dim() << "\n";

        //SPRINKLED COORDINATES
        std::cout << "\nSprinkled values after "<<card<<" sprinklings"<<std::endl;
        std::vector<double> comparisons (6, radius); //for ball and bicone
        if ((strcmp(shape._name, "cylinder")==0))
        {
            comparisons[0] = std::sqrt(radius*radius+myduration/2*myduration/2);
            comparisons[1] = radius;
            comparisons[2] = myduration/2;
            comparisons[3] = myduration/2;
            comparisons[4] = radius;
            comparisons[5] = radius;
        }
        else if ((strcmp(shape._name, "cube")==0))
        {
            comparisons[0] = std::sqrt(dim)*edge/2;
            comparisons[1] = std::sqrt(dim-1)*edge/2;
            comparisons[2] = edge/2;
            comparisons[3] = edge/2;
            comparisons[4] = edge/2;
            comparisons[5] = edge/2;
        }
        else if ((strcmp(shape._name, "cuboid")==0))
        {
            double c0 = edges[0]/2*edges[0]/2; double c1=0;
            for (int i = 0; i<dim; i++)
            {
                double e2 = edges[i]/2*edges[i]/2;
                c0 += e2;
                c1 += e2;
            }
            comparisons[0] = std::sqrt(c0);
            comparisons[1] = std::sqrt(c1);
            comparisons[2] = edges[0]/2;
            comparisons[3] = edges[0]/2;
            comparisons[4] = edges[1]/2;
            comparisons[5] = edges[1]/2;
        }
        std::cout<<"Max Eu Distance: "<<C.max_eu_dist()<<" <= "
                 <<comparisons[0] <<"        ->  "
                 <<(C.max_eu_dist()<=comparisons[0])
                 <<std::endl;
        std::cout<<"Max Sp Radius  : "<<C.max_sp_rad() <<" <= "
                 <<comparisons[1]<<"         ->  "
                 <<(C.max_sp_rad()<=comparisons[1])
                 <<std::endl;
        std::cout<<"Max Time       : "<<C.max_along(0) <<" <= "
                 <<comparisons[2]<<"         ->  "
                 <<(C.max_along(0)<=comparisons[2])
                 <<std::endl;
        std::cout<<"Min Time       : "<<C.min_along(0) <<"  >= -"
                 <<comparisons[3]<<"       ->  "
                 <<(C.min_along(0)<=comparisons[3])
                 <<std::endl;
        std::cout<<"Max Along x    : "<<C.max_along(1) <<" <= "
                 <<comparisons[4]<<"         ->  "
                 <<(C.max_along(1)<=comparisons[4])
                 <<std::endl;
        std::cout<<"Min Along x    : "<<C.min_along(1) <<"  >= -"
                 <<comparisons[5]<<"       ->  "
                 <<(C.min_along(1)<=comparisons[5])
                 <<std::endl;
        // std::cout<<"Max Along y    : "<<C.max_along(2) <<std::endl;
        // std::cout<<"Min Along y    : "<<C.min_along(2) <<std::endl;
        // std::cout<<"Max Along z    : "<<C.max_along(3) <<std::endl;
        // std::cout<<"Min Along z    : "<<C.min_along(3) <<std::endl;

        //CAUSALITY
        auto aretimelike = S.Causality();
        std::vector<double> ovec (dim, 0.0);
        std::vector<double> xvec (dim, 0.0);
        std::vector<double> yvec (dim, 0.0);
        xvec[1] += 1;
        yvec[0] += 1;
        std::cout<<"\nTESTING CAUSALITY\nShould be :          0";
        std::cout<<"\nFrom Spacetime:      "<<aretimelike(ovec, xvec, {}, 0)[0];
        std::cout<<"\nFrom EmbeddedCauset: "<<C.areTimelike(ovec, xvec);
        std::cout<<"\nShould be :          1";
        std::cout<<"\nFrom Spacetime:      "<<aretimelike(ovec, yvec, {}, 0)[0];
        std::cout<<"\nFrom EmbeddedCauset: "<<C.areTimelike(ovec, yvec);
        std::cout<<"\nShould be :          1";
        std::cout<<"\nFrom Spacetime:      "<<aretimelike(yvec, ovec, {}, 0)[0];
        std::cout<<"\nFrom EmbeddedCauset: "<<C.areTimelike(yvec, ovec);

        if (want_coords){
        std::cout << "\nCoordinates:\n";
        print_vector(C._coords);}

        if (want_matrix){
        std::cout<<"\nCausal Matrix:\n";
        print_vector(C._CMatrix);}

        std::cout << "Doing MMd....." << std::endl;
        // MMd estimation
        std::vector<double> MMd_result = C.MMdim_est("big",20,1000);
        std::cout << "MMd (mean,std) = ";
        print_vector(MMd_result);
    }


    std::cout << "\n =========This file works!==========\n" << std::endl;


    auto stop = high_resolution_clock::now();
    double duration = duration_cast<microseconds>(stop - start).count();
    std::cout << "\nTime taken: "
            << duration/pow(10,6) << " seconds" << std::endl;
    return 0;
};