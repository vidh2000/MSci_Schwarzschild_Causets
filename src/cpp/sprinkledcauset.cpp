/// \authors Vid Homsak, Stefano Veroni
/// \date 31/10/2022

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
#include <string.h>
#include <time.h>
#include <vector>
#include <unordered_set>

//using namespace std::chrono;

#include "embeddedcauset.h"
#include "sprinkledcauset.h"
#include "functions.h"
#include "vecfunctions.h"
//#include "causet.h"
//#include "shapes.h"
//#include "spacetimes.h"

using std::vector;
using std::set;
using std::unordered_set;


/**
 * @brief Construct a new Sprinkled Causet object. Note: if shape has 
 * _deltaphi in (0,2pi), then it will restrict to that deltaphi region,
 * where phi is the angle on the xy plane.
 * 
 * @param card: int, number of events to sprinkle
 * @param spacetime: FlatSpacetime.
 * @param shape: CoordinateShape.
 * @param poisson: bool, if true card is used as average of a Poisson
 * distribution to get a new cardinality, otherwise card is used.
 * @param make_matrix: bool, if true make matrix.
 * @param special: bool, if true (and make_matrix) make special Cmatrix
 * @param use_transitivity: bool, if true use also transitivity to establish
 * causality relations. 
 * @param make_sets: bool, if true make set (see sets_type)
 * @param make_links: bool, if true make links (see sets_type)
 * @param sets_type: const char* specifying the type of set:
 * - "past": make _past_links
 * - "future": make _future_links
 * @param seed: int, for coordinate generation. Default 0, i.e. no seed.
 * @param seed: int, for N generation ffrom Poisson. Default 0, i.e. no seed.
 */
SprinkledCauset::SprinkledCauset(int card,
                                 Spacetime spacetime, 
                                 CoordinateShape shape, 
                                 bool poisson,// = true,
                                 bool make_matrix,// = true,
                                 bool special,// = false,
                                 bool use_transitivity,// = true,
                                 bool make_sets,// = false,
                                 bool make_links,// = false,
                                 const char* sets_type,// = "past",
                                 int seed, //=0
                                 int poisson_seed) //=0
                                : EmbeddedCauset()
{
    _spacetime = spacetime;
    _shape = shape; 

    if (poisson)
        {_intensity = card*1;}
    _coords = sprinkle(card, shape, poisson, seed, poisson_seed);

    _size = _coords.size();

    if (strcmp(spacetime._name, "BlackHole")==0){
        Spacetime::CarttoSpherical(_coords);
    }

    this->sort_coords(0, false);
    this->make_attrs("coordinates", make_matrix, special, use_transitivity,
                     make_sets, make_links, sets_type);
}

/**
 * @brief Construct a new Sprinkled Causet object. Note: if shape has 
 * _deltaphi in (0,2pi), then it will restrict to that deltaphi region,
 * where phi is the angle on the xy plane.
 * 
 * @param card: int, number of events to sprinkle
 * @param spacetime: FlatSpacetime.
 * @param shape: CoordinateShape.
 * @param poisson: bool, if true card is used as average of a Poisson
 * distribution to get a new cardinality, otherwise card is used.
 * @param make_matrix: bool, if true make matrix.
 * @param special: bool, if true (and make_matrix) make special Cmatrix
 * @param use_transitivity: bool, if true use also transitivity to establish
 * causality relations. 
 * @param make_sets: bool, if true make set (see sets_type)
 * @param make_links: bool, if true make links (see sets_type)
 * @param sets_type: const char* specifying the type of set:
 * - "past": make _past_links
 * - "future": make _future_links
 * @param seed: int, for coordinate generation. Default 0, i.e. no seed.
 * @param seed: int, for N generation ffrom Poisson. Default 0, i.e. no seed.
 */
SprinkledCauset::SprinkledCauset(double card,
                                 Spacetime spacetime, 
                                 CoordinateShape shape, 
                                 bool poisson,// = true,
                                 bool make_matrix,// = true,
                                 bool special,// = false,
                                 bool use_transitivity,// = true,
                                 bool make_sets,// = false,
                                 bool make_links,// = false,
                                 const char* sets_type,// = "past",
                                 int seed, //=0
                                 int poisson_seed) //=0
                                : EmbeddedCauset()
{
    _spacetime = spacetime;
    _shape = shape; 

    if (poisson)
        {_intensity = card*1;}
    _coords = sprinkle(card, shape, poisson, seed, poisson_seed);

    _size = _coords.size();

    if (strcmp(spacetime._name, "BlackHole")==0){
        Spacetime::CarttoSpherical(_coords);
    }

    this->sort_coords(0, false);
    this->make_attrs("coordinates", make_matrix, special, use_transitivity,
                     make_sets, make_links, sets_type);
}



/**
 * @brief Determines number of points to sprinkle via Poisson distribution or
 *        set as a constant integer and calls the sprinkle_coords method with
 *        the given parameters.
 * 
 * @param card : number (or average) number of points to sprinkle.
 * @param shape : CoordinateShape where to sprinkle
 * @param poisson : bool, if true count is average of Poisson distribution
 * from which get the effective count
 * @param seed : int. Seed for coordinates. Default 0 means no seed.
 * @param poisson_seed : int. Seed for Poisson. Default 0 means no seed.
 * @return vector<vector<double>> coordinates, with ith entry being coordinates
 * of ith point.
 */
vector<vector<double>> SprinkledCauset::sprinkle( int card, 
                                                    CoordinateShape shape,
                                                    bool poisson,// = true,
                                                    int seed, // = 0,
                                                    int poisson_seed)// = 0)
{
    if (card<0)
    {   
        std::cout<<"The sprinkle cardinality has to be a\
            non-negative integer."<<std::endl;
        throw std::invalid_argument(
            "The sprinkle cardinality has to be a non-negative integer.");
    }
    if (poisson)
    {
        if (poisson_seed==0){
            poisson_seed = std::chrono::system_clock::now()
                                                    .time_since_epoch().count();
        }
        
        std::mt19937 gen(poisson_seed);
        std::poisson_distribution<int> distribution(card);
        card = distribution(gen);
    }
    vector<vector<double>> coords = SprinkledCauset::sprinkle_coords
                                                      (card,shape,seed); 
    return coords;
} 


/**
 * @brief Determines number of points to sprinkle via Poisson distribution or
 *        set as a constant integer and calls the sprinkle_coords method with
 *        the given parameters.
 * 
 * @param card : number (or average) number of points to sprinkle.
 * @param shape : CoordinateShape where to sprinkle
 * @param poisson : bool, if true count is average of Poisson distribution
 * from which get the effective count
 * @param seed : int. Seed for coordinates. Default 0 means no seed.
 * @param poisson_seed : int. Seed for Poisson. Default 0 means no seed.
 * @return vector<vector<double>> coordinates, with ith entry being coordinates
 * of ith point.
 */
vector<vector<double>> SprinkledCauset::sprinkle( double card, 
                                                    CoordinateShape shape,
                                                    bool poisson,// = true,
                                                    int seed, // = 0,
                                                    int poisson_seed)// = 0)
{
    if (card<0)
    {   
        std::cout<<"The sprinkle cardinality has to be a\
            non-negative integer."<<std::endl;
        throw std::invalid_argument(
            "The sprinkle cardinality has to be a non-negative integer.");
    }
    if (poisson)
    {
        if (poisson_seed==0){
            poisson_seed = std::chrono::system_clock::now()
                                                    .time_since_epoch().count();
        }
        
        std::mt19937 gen(poisson_seed);
        std::poisson_distribution<int> distribution(card);
        card = distribution(gen);
    }
    vector<vector<double>> coords = SprinkledCauset::sprinkle_coords
                                                      (card,shape,seed); 
    return coords;
}   


/**
 * @brief Generate or ~sprinkle~ "count" coordinates into some "shape"
 * 
 * @param count : fixed number number of points to sprinkle.
 * @param shape : CoordinateShape where to sprinkle
 * @param seed : int, default is no seed.
 * @return vector<vector<double>> coordinates, with ith entry being coordinates
 * of ith point.
 */
vector<vector<double>> SprinkledCauset::sprinkle_coords(int count, 
                                                        CoordinateShape shape,
                                                        int seed)// = 0)
{

    if (count < 0)
    {   
        std::cout<<"The sprinkle cardinality has to be a\
             non-negative integer."<<std::endl;
        throw std::invalid_argument(
            "The sprinkle cardinality has to be a non-negative integer.");
    }

    vector<vector<double>> coords(count,(vector<double>(shape._dim,0.0)));

    // Define mersenne_twister_engine Random Generator (with seed)
    if (seed==0)
    {
        std::random_device rd;
        seed = std::chrono::system_clock::now().time_since_epoch().count(); 
        //seed = rd();
    }  
    std::mt19937 gen(seed);


    if (strcmp(shape._name, "cuboid")==0)
    {
        //vector<double> edges(shape._dim);
        for (int mu = 0; mu<shape._dim; mu++)
        {
            double edge_mu = 0;
            // 1. Get edge value <-- Copied from CoordinateShape::Eges()
            std::string key_mu = "edge_"+std::to_string(mu);
            edge_mu = shape.Parameter(key_mu);
            // 2. // 2. Generate Coordinates (by first finding limits)
            double low  = shape._center[mu] - edge_mu/2.0;
            double high = shape._center[mu] + edge_mu/2.0;
            std::uniform_real_distribution<> dis(low,high);
            for (int n=0; n<count; n++)
            {
                coords[n][mu] = dis(gen);
            }
        }
    }

    else if (strcmp(shape._name, "cube")==0)
    {
        // 1. Get Edge Value
        double edge = shape.Parameter("edge");
        // 2. Generate Coordinates (by first finding limits)
        for (int mu = 0; mu<shape._dim; mu++)
        {
            double low  = shape._center[mu] - edge/2;
            double high = shape._center[mu] + edge/2;
            std::uniform_real_distribution<> dis(low,high);
            for (int n=0; n<count; n++)
            {
                coords[n][mu] = dis(gen);
            }
        }
    }

    //IF IT INVOLVES "CIRCULAR SHAPES"
    else if ((strcmp(shape._name, "ball")==0)     ||
             (strcmp(shape._name, "cylinder")==0) ||
             (strcmp(shape._name, "diamond")==0)  ||
             (strcmp(shape._name, "bicone")==0))
    {
        bool isCylindrical = strcmp(shape._name, "cylinder")==0;
        bool isDiamond     = strcmp(shape._name, "diamond") == 0 ||
                             strcmp(shape._name, "bicone") == 0;
        
        double d = shape._dim;
        // Get radius. Other methods failed.. loop is supershort and unrepeatd
        double R = shape.Parameter("radius");

        if (d==2 && isDiamond)
        {
            //pick "count" random coordinate tuples uniformly:
            vector<vector<double>> uv;
            std::uniform_real_distribution<> dis(-1.0,1.0);
            for (int n = 0; n<count; n++)
            {
                double u = dis(gen);
                double v = dis(gen);
                coords[n][0] = (u + v)*R/2.0;
                coords[n][1] = (u - v)*R/2.0;
            } 
        }

        else
        {
            // number of "different" and "radius-like" coordinates:
            // - 0 for ball
            // - 1 (time) for others
            double n_different = (strcmp(shape._name, "ball")==0)? 0 : 1;
            double n_rad = d - n_different;

            //Get lower bound for hollow case
            double r_low = std::pow(shape.Parameter("hollow"), n_rad);

            //Get time done separately for cylinder
            if (isCylindrical)
            {   
                double duration = shape.Parameter("duration");
                double t_low = shape._center[0]-duration/2;
                double t_high= shape._center[0]+duration/2;
                std::uniform_real_distribution<> cyltime(t_low, t_high);
                for (int n = 0; n<count; n++)
                    {coords[n][0] = cyltime(gen);} 
            }

            // Prepare scaling factor, where to accumulate all
            // scaling factors to apply at the end
            double r_scaling_i = 0;
            std::uniform_real_distribution<> r_uni(r_low, 1);
            std::uniform_real_distribution<> other_uni(0, 1);
            for (int i = 0; i<count; i++)
            {
                //1)
                //First generate random numbers on (n_rad-1)-sphere of rad 1
                std::normal_distribution<> gauss(0,1);
                double r2_gauss = 0;
                for (int a = n_different; a<d; a++)
                {
                    double x_a = gauss(gen);
                    coords[i][a] = x_a;
                    r2_gauss += x_a*x_a;
                }
                // 1/r2 brings the points on a sphere of radius 1
                r_scaling_i = 1/std::sqrt(r2_gauss);
                
                // 2)
                // Scale points radially inside n_rad-ball
                r_scaling_i *= R * std::pow(r_uni(gen), 1/n_rad);
                if (isDiamond)
                {
                    //Set time properly
                    double droot_k = std::pow(other_uni(gen), 1/d);
                    double t_sign = (other_uni(gen)<0.5)? -1. : 1.;
                    coords[i][0] = t_sign * (1-droot_k) * R + shape._center[0];
                    r_scaling_i *= droot_k;
                }
                for (int a = n_different; a<d; a++)
                {
                    coords[i][a] *= r_scaling_i;
                    coords[i][a] += shape._center[a];
                }

                // If deltaphi in (0, 2pi)
                // place points inside slice of said deltaphi
                double deltaphi = shape._deltaphi;
                if (deltaphi > 0 && deltaphi <= 2*3.1415926535)
                {
                    double phi0 = std::atan2(coords[i][2],coords[i][1]);
                    phi0 += (phi0>0)? 0 : 2*3.1415926535; 
                    double phi = phi0 * deltaphi / (2*3.1415926535);

                    double rho = std::sqrt(coords[i][1]*coords[i][1]
                                          +coords[i][2]*coords[i][2]);
                    coords[i][1] = rho * std::cos(phi);
                    coords[i][2] = rho * std::sin(phi);
                }
            }
        }
    }
    return coords;
}



// Run with
// cd scripts_cpp/causet_cpp
// g++ -g causet.cpp shapes.cpp spacetimes.cpp sprinkledcauset.cpp embeddedcauset.cpp -std=c++17 -o sprinkledcauset -O2
// .\sprinklededcauset
// rm sprinkledcauset.exe
// cd ../
// cd../
// int main(){
//     std::cout << "sprinkledcauset.cpp WORKS! :)";
// }


