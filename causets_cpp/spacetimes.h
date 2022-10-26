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
        std::map<char, float> _params;

        
        int&  Dim        = _dim;
        char& Name       = _name;
        char& MetricName = _metricname; 
        virtual float Parameter (char key);


        virtual CoordinateShape DefaultShape();  

        virtual std::vector<float>_T_slice_sampling(float t, 
                                                    std::vector<float>origin,
                                                    int samplingsize = -1);     
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
        
        FlatSpacetime(int dim = 2, std::vector<float>period = {0});
        float ds2    (std::vector<float> xvec, std::vector<float> yvec);
        float ds     (std::vector<float> xvec, std::vector<float> yvec);

        float           Parameter (char key);
        CoordinateShape DefaultShape();
        
        typedef std::vector<bool> (*func)
        (std::vector<float> xvec, std::vector<float> yvec);
        func Causality();       
};


/*=============================================================================
=============================================================================*/
class BlackHoleSpacetime: public Spacetime
/**
 * @brief Implementation of black hole spacetimes, globally hyperbolic.
 */
{
    public:
        
        BlackHoleSpacetime(int dim = 2, float r_S = 0.5,
                            char metric = 'Eddington-Finkelstein');
        float ds2    (std::vector<float> xvec, std::vector<float> yvec);
        float ds     (std::vector<float> xvec, std::vector<float> yvec);

        float           Parameter (char key);
        CoordinateShape DefaultShape();
        
        typedef std::vector<bool> (*func)
        (std::vector<float> xvec, std::vector<float> yvec);
        func Causality();     
};
#endif /* SPACETIMES_H */