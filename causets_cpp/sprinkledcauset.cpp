/// \authors Vid Homsak, Stefano Veroni
/// \date 31/10/2022

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <set>
#include <stack>
#include <stdio.h>
#include <stdexcept>
#include <string>
#include <vector>
#include <chrono>
#include <unordered_set>

//using namespace std::chrono;


#include "functions.h"
#include "MyVecFunctions.h"
#include "causet.h"
#include "embeddedcauset.h"
#include "sprinkledcauset.h"
#include "shapes.h"
#include "spacetimes.h"

using std::vector;
using std::set;
using std::unordered_set;

// For Vid (can't locate headers for some reason)
// Path: D:\Documents\Sola\Imperial College London\Year 4\MSci project\
            Project\causets_code\causets_cpp\"header".h...




//////////////////////////////////////////////////////////////////////////////
// ::SprinkledCauset::
// Handles a causal set that is embedded in a subset of a manifold.
//////////////////////////////////////////////////////////////////////////////

/**
 * @brief Construct a new Sprinkled Causet object. 
 * 
 * @param card: int, number of events to sprinkle
 * @param spacetime: Spacetime.
 * @param shape: CoordinateShape.
 * @param poisson: bool, if true card is used as average of a Poisson
 * distribution to get a new cardinality, otherwise card is used.
 * @param make_matrix: bool, if true make matrix.
 * @param special: bool, if true (and make_matrix) make special Cmatrix
 * @param use_transitivity: bool, if true use also transitivity to establish
 * causality relations. 
 * @param make_sets: bool, if true make set (see sets_type)
 * @param make_links: bool, if true make links (see sets_type)
 * @param method: const char* specifying
 * - "intensity": card is the average of Poisson distribution used to extract 
 * number of events to sprinkle
 * - "card": the number of events to sprinkle is fixed to card
 * @param sets_type: const char* specifying the type of set:
 * - "past": make _past_links
 * - "future": make _future_links
 * @param seed: int, for coordinate generation.
 */
SprinkledCauset::SprinkledCauset(int card,
                                 Spacetime spacetime, 
                                 CoordinateShape shape, 
                                 bool poisson = false,
                                 bool make_matrix = true,
                                 bool special = false,
                                 bool use_transitivity = true,
                                 bool make_sets = false,
                                 bool make_links = false,
                                 const char* sets_type = "past",
                                 int seed = std::nanf(""))
{
    if (poisson)
        {_intensity = card*1;}
    
    _coords = sprinkle_coords(card, shape, seed);
    _size = _coords.size();
    _spacetime = spacetime;
    _shape = shape; 
    
    if (make_matrix)
        {this->make_cmatrix("coordinates", special, use_transitivity,
                                make_sets, make_links, sets_type);}
    
    else if (make_sets)
        {
            if (make_links)
            {
                if (sets_type == "past")
                    {this->make_all_pasts("coordinates");}
                else if (sets_type == "future")
                    {this->make_all_futures("coordinates");}
            }
            else
            {
                if (sets_type == "past")
                    {this->make_pasts("coordinates");}
                else if (sets_type == "future")
                    {this->make_futures("coordinates");}
            }
        }

    else if (make_links)
        {
            if (sets_type == "past")
                    {this->make_past_links("coordinates");}
            else if (sets_type == "future")
                {this->make_fut_links("coordinates");}   
        }
        
    else
        {
            throw std::invalid_argument("At least one among make_matrix, \
                                    make_sets and make_links must be true");
        }
}


// Methods
static vector<vector<double>> SprinkledCauset::sprinkle_coords(
                                            int count, CoordinateShape shape,
                                            int seed = std::nanf(""))
{
    if (count < 0)
    {   throw std::invalid_argument(
            "The sprinkle cardinality has to be a non-negative integer.");
    }
    vector<vector<double>> coords;
    if ((shape._name == "cube") || (shape._name == "cuboid"))
    {
        vector<double> low;
        vector<double> high;
    
        low = shape.Center()-shape.Edges()/2;
        high = shape.Center()-shape.Edges()/2;

        // Generate coords randomly

        // Will be used to obtain a seed for the random number engine
        if (std::isnan(seed))
            {
                std::random_device rd;
                seed = rd();
            }  
        // Standard mersenne_twister_engine seeded with rd()
        std::mt19937 gen(seed); 
        std::uniform_real_distribution<> dis(low,high);
        for (int i=0; i<count;i++)
        {
            coords[i,:] = dis(gen);
        }
    }
    else if ((shape._name == "ball") || (shape._name == "cylinder") ||
             (shape._name == "diamond") || (shape._name == "bicone"))
    {
        // Create circle based sprinkle:
        bool isCylindrical = shape.Name()=="cylinder";
        bool isDiamond = (shape.Name()=="diamond") || (shape.Name()=="bicone")

        int d = *this.Dim()
        double b_r = shape.Parameter().radius;
        if (d==2 && isDiamond)
        {
            //pick "count" random coordinate tuples uniformly:
            vector<vector<double>> uv;
            // Random generator stuff
            std::random_device rd;  
            std::mt19937 gen(rd()); 
            std::uniform_real_distribution<> dis(-1.0,1.0);
            for (i=0;i<count;i++)
            {
                vector<double> = 
                for (j=0;j<2;j++)
                {
                    
                    /* placeholder. Find out how to generate a 
                    random matrix (count,2) size in c++*/ 
                }
                
            } 
        }
    }
    
}

static vector<vector<double>> SprinkledCauset::sprinkle( int count, 
                                                    CoordinateShape shape,
                                                    bool poisson = false,
                                                    int seed = std::nanf(""))
{
    if (count<0)
        {   throw std::invalid_argument(
                "The sprinkle cardinality has to be a non-negative integer.");
        }
    if (poisson)
        {
            std::default_random_engine generator;
            std::poisson_distribution<int> distribution(count);
            count = distribution(generator);
        }
    vector<vector<double>> coords = SprinkledCauset::sprinkle_coords
                                                      (count,shape,seed); 
    return coords;
}   





