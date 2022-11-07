/// \authors Vid Homsak, Stefano Veroni
/// \date 25/09/2022

#ifndef SPACETIMES_H
#define SPACETIMES_H

#include <algorithm>
#include <cmath>
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
        const char* _metricname;

        bool _isPeriodic = false;
        std::vector<double> _period = {};

        // CONSTRUCTOR
        Spacetime();

        // VIRTUAL FUNCTIONS
      
        std::vector<double>_T_slice_sampling(double t, 
                                                    std::vector<double>origin,
                                                    int samplingsize = -1); 

        // CAUSALITY
        typedef std::vector<bool> (*func)
        (std::vector<double> xvec, std::vector<double> yvec, 
         std::vector<double> period);
        virtual func Causality();  
        
        static std::vector<bool> causal1d(std::vector<double> xvec, 
                                          std::vector<double> yvec,
                                          std::vector<double> period={});
        
        ~Spacetime();
};



/*=============================================================================
=============================================================================*/
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
        
        FlatSpacetime(int dim = 4, std::vector<double> period = {});

        double ds2    (std::vector<double> xvec, std::vector<double> yvec);
        double ds     (std::vector<double> xvec, std::vector<double> yvec);
        
        typedef std::vector<bool> (*func)
        (std::vector<double> xvec, std::vector<double> yvec, std::vector<double> period);
        func Causality();   

        static std::vector<bool> causal (std::vector<double> xvec, 
                                    std::vector<double> yvec,
                                    std::vector<double> period={});
        static std::vector<bool> causal_periodic (std::vector<double> xvec, 
                                            std::vector<double> yvec,
                                            std::vector<double>period={});
        
        ~FlatSpacetime();
};


/*=============================================================================
=============================================================================*/
class BlackHoleSpacetime: public Spacetime
/**
 * @brief Implementation of black hole spacetimes, globally hyperbolic.
 */
{
    public:
        
        // Constructor
        BlackHoleSpacetime(int dim = 2,
                           double r_S = 0.5,
                           const char* metric = "Eddington-Finkelstein");
                
        // Methods
        double ds2(std::vector<double> xvec, std::vector<double> yvec);

        double ds(std::vector<double> xvec, std::vector<double> yvec);
        
        typedef std::vector<bool> (*func)
        (std::vector<double> xvec, std::vector<double> yvec, 
        std::vector<double> period);
        func Causality();   

        ~BlackHoleSpacetime();  
};
#endif /* SPACETIMES_H */