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


// Sprinkled causet parameters
int card = 200;
int dim = 3;
std::vector<double> center (dim, 0.0);
double radius = 1.5;
double myduration = 5;
double hollow = 0.7;
double edge = 1.5;
std::vector<double> edges = {1,2,3,4};

// For BH
double mass = 0.5;


// Sprinkling Parameters
bool poisson = false;
bool make_matrix = true;
bool special = false;
bool use_transitivity = false;
bool make_sets = false;
bool make_links = false; 
const char* sets_type = "future"; 
const char* name = "cylinder";

//Analysis Parameters (printing related)
bool want_matrix = false;
bool want_coords = false;


int main(){
    
    auto start = high_resolution_clock::now();
    //std::cout << "Starting building shape..." << std::endl;
    std::vector<const char*> names = {"cylinder"};// "ball", "cylinder", "cube", "cuboid"};
    edges.resize(dim);
    //card = (want_coords || want_matrix)? 10 : card;
    for (const char* name : names)
    {

        CoordinateShape shape(dim,name,center,radius,myduration,hollow,edge,edges);
        std::cout<<"\n\n============= USING "<<name<<" ====================\n";
        // std::cout << "What are the just-Shape's parameters?\n";
        // for (auto const& p : shape._params)
        //     {std::cout << "-- " << p.first << "=" << p.second << '\n';}
        

        Spacetime S = Spacetime();
        //S.FlatSpacetime(dim);
        S.FlatSpacetime(dim);
        SprinkledCauset C(card, S, shape, poisson,
                          make_matrix, special, use_transitivity,
                          make_sets, make_links,sets_type);

        //PARAMETERS
        //std::cout << "\nWhat are the Causet's Shape's parameters at the end?\n";
        //for (auto const& p : C._shape._params){
        //    std::cout << "-- " << p.first << '=' << p.second << '\n';}
        // if ((strcmp(shape._name, "ball")==0)      ||
        //      (strcmp(shape._name, "cylinder")==0) ||
        //      (strcmp(shape._name, "diamond")==0)  ||
        //      (strcmp(shape._name, "bicone")==0))
        // {
        //     std::cout<<"Radius obtained from Parameter: "
        //          <<C._shape.Parameter("radius")<<std::endl;
        // }
        // std::cout<<"-- cardinality=" << C._size << "\n";
        // std::cout<<"-- dim=" << C.spacetime_dim() << "\n";

        //SPRINKLED COORDINATES
        // std::cout << "\nSprinkled values after "<<card<<" sprinklings"<<std::endl;
        // std::vector<double> comparisons (6, radius); //for ball and bicone
        // if ((strcmp(shape._name, "cylinder")==0))
        // {
        //     comparisons[0] = std::sqrt(radius*radius+myduration/2*myduration/2);
        //     comparisons[1] = radius;
        //     comparisons[2] = myduration/2;
        //     comparisons[3] = myduration/2;
        //     comparisons[4] = radius;
        //     comparisons[5] = radius;
        // }
        // else if ((strcmp(shape._name, "cube")==0))
        // {
        //     comparisons[0] = std::sqrt(dim)*edge/2;
        //     comparisons[1] = std::sqrt(dim-1)*edge/2;
        //     comparisons[2] = edge/2;
        //     comparisons[3] = edge/2;
        //     comparisons[4] = edge/2;
        //     comparisons[5] = edge/2;
        // }
        // else if ((strcmp(shape._name, "cuboid")==0))
        // {
        //     double c0 = edges[0]/2*edges[0]/2; double c1=0;
        //     for (int i = 0; i<dim; i++)
        //     {
        //         double e2 = edges[i]/2*edges[i]/2;
        //         c0 += e2;
        //         c1 += e2;
        //     }
        //     comparisons[0] = std::sqrt(c0);
        //     comparisons[1] = std::sqrt(c1);
        //     comparisons[2] = edges[0]/2;
        //     comparisons[3] = edges[0]/2;
        //     comparisons[4] = edges[1]/2;
        //     comparisons[5] = edges[1]/2;
        // }
        // std::cout<<"Max Eu Distance: "<<C.max_eu_dist()<<" <= "
        //          <<comparisons[0] <<"        ->  "
        //          <<(C.max_eu_dist()<=comparisons[0])
        //          <<std::endl;
        // std::cout<<"Max Sp Radius  : "<<C.max_sp_rad() <<" <= "
        //          <<comparisons[1]<<"         ->  "
        //          <<(C.max_sp_rad()<=comparisons[1])
        //          <<std::endl;
        // std::cout<<"Max Time       : "<<C.max_along(0) <<" <= "
        //          <<comparisons[2]<<"         ->  "
        //          <<(C.max_along(0)<=comparisons[2])
        //          <<std::endl;
        // std::cout<<"Min Time       : "<<C.min_along(0) <<"  >= -"
        //          <<comparisons[3]<<"       ->  "
        //          <<(C.min_along(0)<=comparisons[3])
        //          <<std::endl;
        // std::cout<<"Max Along x    : "<<C.max_along(1) <<" <= "
        //          <<comparisons[4]<<"         ->  "
        //          <<(C.max_along(1)<=comparisons[4])
        //          <<std::endl;
        // std::cout<<"Min Along x    : "<<C.min_along(1) <<"  >= -"
        //          <<comparisons[5]<<"       ->  "
        //          <<(C.min_along(1)<=comparisons[5])
        //          <<std::endl;


        // std::cout<<"Max Along y    : "<<C.max_along(2) <<std::endl;
        // std::cout<<"Min Along y    : "<<C.min_along(2) <<std::endl;
        // std::cout<<"Max Along z    : "<<C.max_along(3) <<std::endl;
        // std::cout<<"Min Along z    : "<<C.min_along(3) <<std::endl;

        // //CAUSALITY
        // auto aretimelike = S.Causality();
        // std::vector<double> ovec (dim, 0.0);
        // std::vector<double> xvec (dim, 0.0);
        // std::vector<double> yvec (dim, 0.0);
        // xvec[1] += 1;
        // yvec[0] += 1;
        // std::cout<<"\nTESTING CAUSALITY\nShould be :          0";
        // std::cout<<"\nFrom Spacetime:      "<<aretimelike(ovec, xvec, {}, 0)[0];
        // std::cout<<"\nFrom EmbeddedCauset: "<<C.areTimelike(ovec, xvec);
        // std::cout<<"\nShould be :          1";
        // std::cout<<"\nFrom Spacetime:      "<<aretimelike(ovec, yvec, {}, 0)[0];
        // std::cout<<"\nFrom EmbeddedCauset: "<<C.areTimelike(ovec, yvec);
        // std::cout<<"\nShould be :          1";
        // std::cout<<"\nFrom Spacetime:      "<<aretimelike(yvec, ovec, {}, 0)[0];
        // std::cout<<"\nFrom EmbeddedCauset: "<<C.areTimelike(yvec, ovec);

        if (want_coords){
        std::cout << "\nCoordinates:\n";
        print_vector(C._coords);}

        if (want_matrix){
        std::cout<<"\nCausal Matrix:\n";
        print_vector(C._CMatrix);}

        // std::cout << "C_pasts, N_pasts = " << C._pasts.size() << std::endl;
        // for (auto past : C._pasts)
        // {
        //     print_set(past);
        // }
        // std::cout << "C_futures, N_futures = " << C._futures.size() << std::endl;
        // for (auto fut : C._futures)
        // {
        //     print_set(fut);
        // }

        // std::cout << "\nDoing MMd....." << std::endl;
        // // MMd estimation
        // std::vector<double> MMd_result = C.MMdim_est("big",20,
        //                             vecmin(std::vector<int> {1000,C._size/2}));
        // std::cout << "MMd (mean,std) = ";
        // print_vector(MMd_result);
    }


    std::cout << "\n =========This file works!==========\n" << std::endl;


    auto stop = high_resolution_clock::now();
    double duration = duration_cast<microseconds>(stop - start).count();
    std::cout << "Time taken for N= " << card << ", dim = "<< dim << ": "
            << duration/pow(10,6) << " seconds" << std::endl;
    
    return 0;
};