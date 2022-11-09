/// \authors Vid Homsak, Stefano Veroni
/// \date 25/09/2022

#ifndef SPACETIMES_H
#define SPACETIMES_H

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>
#include <fstream>
#include <stack>
#include <string>
#include <stdio.h>
#include <vector>

#include "shapes.h"


class Spacetime
/**
 * @brief Super class for implementation of spacetimes
 */
{
    public:
        int  _dim;
        const char* _name;
        std::string _metricname;
        double _mass = 0; //for Black Hole, but required here for generality

        bool _isPeriodic = false;
        std::vector<double> _period = {};

        // CONSTRUCTOR
        Spacetime();

        // FUNCTIONS
        std::vector<double>_T_slice_sampling(double t, 
                                                    std::vector<double>origin,
                                                    int samplingsize = -1); 

        // CAUSALITY
        typedef std::vector<bool> (*func)
        (std::vector<double> xvec, std::vector<double> yvec, 
        std::vector<double> period, double mass);
        virtual func Causality();  
        
        static std::vector<bool> causal1d(std::vector<double> xvec, 
                                          std::vector<double> yvec,
                                          std::vector<double> period={},
                                          double mass = 0);
        
        //virtual ~Spacetime();
};



/*=============================================================================
=============================================================================*/


/**
 * @brief Construct a new Flat Spacetime:: Flat Spacetime object
 * 
 * @param dim int:   dimension of flat spacetime. Default 4.
 * @param period vector<double>:   periodicity vector: ith entry is period 
 * along ith SPATIAL diemnsion (0->x, 1->y, 2->z). Default is no periodicty.
 * No periodicity along ith direction requires 0. Usually, but not necessarily,
 * with cuboid might want period[i]=shape_object.Edges()[i].
 */
class FlatSpacetime: public Spacetime
/**
 * @brief Initializes Minkowski spacetime for dim >= 1.
    As additional parameter, the spatial periodicity can be specified (using 
    the key "period") as float (to be applied for all spatial directions 
    equally) or as tuple (with a float for each spatial dimension). A positive 
    float switches on the periodicity along the respective spatial direction, 
    using this value as period. The default is 0.0, no periodicity in any 
    direction.
 */
{
    public:
        // int  _dim;
        // const char* _name;
        // std::string _metricname;
        
        FlatSpacetime(int dim = 4, std::vector<double> period = {});

        double ds2    (std::vector<double> xvec, std::vector<double> yvec);
        double ds     (std::vector<double> xvec, std::vector<double> yvec);
        
        typedef std::vector<bool> (*func)
        (std::vector<double> xvec, std::vector<double> yvec, 
        std::vector<double> period, double mass);
        func Causality();   

        static std::vector<bool> causal (std::vector<double> xvec, 
                                            std::vector<double> yvec,
                                            std::vector<double> period={}, 
                                            double mass = 0);
        static std::vector<bool> causal_periodic (std::vector<double> xvec, 
                                            std::vector<double> yvec,
                                            std::vector<double>period={}, 
                                            double mass = 0);
        
        //~FlatSpacetime();
};


/*=============================================================================
=============================================================================*/
class BlackHoleSpacetime: public Spacetime
/**
 * @brief Implementation of black hole spacetimes, globally hyperbolic.
 */
{
    public:
        // int  _dim;
        // const char* _name;
        // std::string _metricname;
        double _r_S;

        // Constructor
        BlackHoleSpacetime(int dim = 2,
                           double r_S = 0.5,
                           std::string metric = "Eddington-Finkelstein");
    
        // Methods
        double ds2(std::vector<double> xvec, std::vector<double> yvec);

        double ds(std::vector<double> xvec, std::vector<double> yvec);
        
        typedef std::vector<bool> (*func)
        (std::vector<double> xvec, std::vector<double> yvec, 
        std::vector<double> period, double mass);
        func Causality();   

        static std::vector<bool> causal (std::vector<double> xvec, 
                                         std::vector<double> yvec,
                                         std::vector<double> period={},
                                         double mass = 0);

        ~BlackHoleSpacetime();  
};


#endif /* SPACETIMES_H */