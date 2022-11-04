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
        char _name;
        char _metricname;
        std::map<char, double> _params;

        //Getters
        int&  Dim        = _dim;
        char& Name       = _name;
        char& MetricName = _metricname; 
        virtual double Parameter (char key);
        virtual CoordinateShape DefaultShape();  
        
        virtual std::vector<double>_T_slice_sampling(double t, 
                                                    std::vector<double>origin,
                                                    int samplingsize = -1); 

        typedef std::vector<bool> (*func)
        (std::vector<double> xvec, std::vector<double> yvec);
        func Causality();  
        
        static std::vector<bool> causal1d(std::vector<double> xvec, 
                                          std::vector<double> yvec);
};



/*=============================================================================
=============================================================================*/
class FlatSpacetime: public Spacetime
/**
 * @brief Initializes Minkowski spacetime for dim >= 1.
    As additional parameter, the spatial periodicity can be specified (using 
    the key 'period') as float (to be applied for all spatial directions 
    equally) or as tuple (with a float for each spatial dimension). A positive 
    float switches on the periodicity along the respective spatial direction, 
    using this value as period. The default is 0.0, no periodicity in any 
    direction.
 */
{
    public:
        bool _isPeriodic;
        
        FlatSpacetime(int dim = 2);
        
        CoordinateShape DefaultShape();

        double ds2    (std::vector<double> xvec, std::vector<double> yvec);
        double ds     (std::vector<double> xvec, std::vector<double> yvec);
        
        typedef std::vector<bool> (*func)
        (std::vector<double> xvec, std::vector<double> yvec);
        func Causality();   

        static std::vector<bool> causal (std::vector<double> xvec, 
                                         std::vector<double> yvec);
};


/*=============================================================================
=============================================================================*/
class BlackHoleSpacetime: public Spacetime
/**
 * @brief Implementation of black hole spacetimes, globally hyperbolic.
 */
{
    public:
        
        BlackHoleSpacetime(int dim = 2, double r_S = 0.5,
                           char metric = 'Eddington-Finkelstein');
        
        
        double ds2    (std::vector<double> xvec, std::vector<double> yvec);

        
        double ds     (std::vector<double> xvec, std::vector<double> yvec);

        CoordinateShape DefaultShape();
        
        typedef std::vector<bool> (*func)
        (std::vector<double> xvec, std::vector<double> yvec);
        func Causality();     
};
#endif /* SPACETIMES_H */