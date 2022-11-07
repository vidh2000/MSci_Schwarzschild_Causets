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

using std::vector;


class Spacetime
/**
 * @brief Super class for implementation of spacetimes
 */
{
    public:
        int  _dim;
        const char* _name;
        const char* _metricname;
        std::map<const char*, double> _params;

        // CONSTRUCTOR
        Spacetime();

        // VIRTUAL GETTERS (BECOME SPECIFIC IN SUBCLASSES)
        virtual double Parameter (const char* key);
        virtual CoordinateShape DefaultShape();  
        
        virtual vector<double>_T_slice_sampling(double t, 
                                                    vector<double>origin,
                                                    int samplingsize = -1); 

        // CAUSALITY
        typedef vector<bool> (*func)
        (vector<double> xvec, vector<double> yvec, vector<double> period);
        func Causality();  
        
        static vector<bool> causal1d(vector<double> xvec, vector<double> yvec, 
                                vector<double>period={});
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
        bool _isPeriodic;
        vector<double> _period = {};
        
        FlatSpacetime(int dim = 4, vector<double> period = {});
        
        CoordinateShape DefaultShape();

        double ds2    (vector<double> xvec, vector<double> yvec);
        double ds     (vector<double> xvec, vector<double> yvec);
        
        typedef vector<bool> (*func)
        (vector<double> xvec, vector<double> yvec, vector<double> period);
        func Causality();   

        static vector<bool> causal (vector<double> xvec, 
                                    vector<double> yvec,
                                    vector<double> period={});
        static vector<bool> causal_periodic (vector<double> xvec, 
                                            vector<double> yvec,
                                            vector<double>period);
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
        double ds2(vector<double> xvec, vector<double> yvec);

        double ds(vector<double> xvec, vector<double> yvec);

        CoordinateShape DefaultShape();
        
        typedef vector<bool> (*func)
        (vector<double> xvec, vector<double> yvec);
        func Causality();     
};
#endif /* SPACETIMES_H */