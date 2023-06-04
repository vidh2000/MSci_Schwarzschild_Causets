#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <stdio.h>
#include <string.h>
#include <stack>
#include <string>
#include <vector>

#include "boost/numeric/odeint.hpp"

#include "shapes.h"
#include "spacetimes.h"
#include "vecfunctions.h"
#include "functions.h"

using std::vector;

#define M_PI  3.14159265358979323846  /* pi */



/*============================================================================
* ============================================================================
* GENERAL SPACETIME METHODS
* ============================================================================
============================================================================*/

/**
 * @brief Construct a new Spacetime:: Spacetime object. It does nothing 
 * but initialise a Spacetime onject. Then call:
 * - FlatSpacetime(dim, period)
 * - BlackHoleSpacetime(dim)
 */
Spacetime::Spacetime(){}


/**
 * @brief Internal function for the time sampling array for a cone from 
 * "origin point"(hence take t0=origin[0]) to time "t".
 */
vector<double> Spacetime::T_slice_sampling(double t, 
                                            vector<double>origin,
                                            int samplingsize)
{
    vector<double> ts (samplingsize);
    double t_i = (t-origin[0])/samplingsize;
    for (int i = 0; i<samplingsize; i++)
        {ts[i] += origin[0] + t_i*i;}
    return ts;
}     



typedef bool (*func)
(const vector<double>& xvec, const vector<double>& yvec, 
 vector<double> period, double mass);
/**
 * @brief Return callable "bool callable (const vector<double>& xvec, 
 * const vector<double>& yvec, vector<double> period, double mass", i.e.
 * returns a function that returns the bool saying if x-y is timelike
 * based on spacetime features.
 */
func Spacetime::Causality()
{
    if (_dim == 1)
        {return &Spacetime::causal1d;}
    else if (std::strcmp(_name, "Flat")==0)
    {
        if (_isPeriodic) return &Spacetime::Flat_causal_periodic;
        else return &Spacetime::Flat_causal;
    }
    else /*(std::strcmp(_name, "BlackHole")==0)*/
    {
        // Changing coordinates into EF (original)
        // place holder //
        if      (_dim == 4) return &Spacetime::BH_causal4D;
        else if (_dim == 3) return &Spacetime::BH_causal3D;
        else if (_dim == 2) return &Spacetime::BH_causal2D;
        else
        {
            std::cout << "Please use D=2, D=3 or D=4 for Schwarzschild\n";
            throw std::invalid_argument("problem");
        }
    }
}


typedef vector<bool> (*general_func)
(const vector<double>& xvec, const vector<double>& yvec, vector<double> period, double mass);
/**
 * @brief Return callable "vector<bool> callable (const vector<double>& xvec, 
 * const vector<double>& yvec, vector<double> period, double mass", i.e.
 * returns a function that returns a vector of bools {x-y timelike, x<=y, x>y},
 * based on spacetime features.
 */
general_func Spacetime::General_Causality()
{
    if (_dim == 1)
        {return &Spacetime::general_causal1d;}
    else if (std::strcmp(_name, "Flat")==0)
    {
        if (_isPeriodic) return &Spacetime::Flat_general_causal_periodic;
        else return &Spacetime::Flat_general_causal;
    }
    else /*(std::strcmp(_name, "BlackHole")==0)*/
    {
        // Changing coordinates into EF (original)
        // place holder //
        if      (_dim == 4) return &Spacetime::BH_general_causal4D;
        else if (_dim == 3) return &Spacetime::BH_general_causal3D;
        else if (_dim == 2) return &Spacetime::BH_general_causal2D;
        else
        {
            std::cout << "Please use D=2, D=3 or D=4 for Schwarzschild\n";
            throw std::invalid_argument("problem");
        }
    }
}


/**
 * @brief Simplest function of two events in 1D
 * 
 * @param xvec : vector<double>(1), time-coordinate of x
 * @param yvec : vector<double>(1), time-coordinate of y
 * @param period : needed for consistency with Causality, not used
 * @return bool : 1
 */
bool Spacetime::causal1d(const vector<double>& xvec, const vector<double>& yvec,
                         vector<double> period, double mass)
{
    return true;
}


/**
 * @brief Simplest function of two events in 1D
 * 
 * @param xvec : vector<double>(1), time-coordinate of x
 * @param yvec : vector<double>(1), time-coordinate of y
 * @param period : needed for consistency with Causality, not used
 * @return vector<bool> : {1, x<=y, x>y}
 */
vector<bool> Spacetime::general_causal1d(const vector<double>& xvec, 
                                         const vector<double>& yvec,
                                         vector<double> period, double mass)
{
    //std::cout<<"For debuggin: calling Spacetime::causal1d\n";
    double t_delta = yvec[0] - xvec[0];
    return {1, t_delta >= 0.0, t_delta < 0.0};
}





/*============================================================================
* ============================================================================
* FLAT SPACETIME
* ============================================================================
============================================================================*/

/**
 * @brief Construct a new Flat Spacetime:: Flat Spacetime object
 * 
 * @param dim int:   dimension of flat spacetime. Default 4.
 * @param period vector<double>:   periodicity vector: ith entry is period 
 * along ith SPATIAL diemnsion (0->x, 1->y, 2->z). Default is no periodicty.
 * No periodicity along ith direction requires 0. Usually, but not necessarily,
 * with cuboid might want period[i]=shape_object.Edges()[i].
 */
void Spacetime::FlatSpacetime(int dim, vector<double> period)
{
    if (dim < 1)
    {
        std::cout<<"Given dim was smaller than 1."<<std::endl;
        throw std::invalid_argument("Dimension has to be at least 1.");
    }
    else if (period.size()!=0 && period.size()!=dim-1)
    {
        std::cout<<"Period size was "<<period.size()<<", it has to be either\
 0 or dim-1. Note: "<<"dimension was "<<dim<<std::endl;
        throw std::invalid_argument("Period's size is wrong!");
    }

    _dim = dim;
    _name = "Flat";
    _metricname = "Minkowski";
    if (period.size() && std::count(period.begin(), period.end(), 0)!=dim-1) 
    {
        _isPeriodic = true;
        _period = period;
        for (int i=0; i<period.size(); i++)
        {
            if (period[i]<0) {_period[i]=-period[i];}
        }
    }
}


/**
 * @brief Spacetime Interval Squared t squared - x_i*x^i
 */
double Spacetime::Flat_ds2(const vector<double>& xvec, const vector<double>& yvec)
{
    
    double time_delta2  = (yvec[0] - xvec[0])*(yvec[0] - xvec[0]);
    double space_delta2 = 0;
    for (int i = 1; i<_dim; i++)
        {space_delta2 += (yvec[i]-xvec[i])*(yvec[i]-xvec[i]);}
    return time_delta2 - space_delta2;
}


/**
 * @brief Spacetime Interval sqrt(t squared - x_i*x^i)
 */
double Spacetime::Flat_ds(const vector<double>& xvec, const vector<double>& yvec)
{
    
    double time_delta2  = (yvec[0] - xvec[0])*(yvec[0] - xvec[0]);
    double space_delta2 = 0;
    for (int i = 1; i<_dim; i++)
        {space_delta2 += (yvec[i]-xvec[i])*(yvec[i]-xvec[i]);}
    return std::sqrt(time_delta2 - space_delta2);
}


/**
 * @brief Function of two events in any D returning true if x-y timelike?,
 * else false.
 * 
 * @param xvec vector<double>:  coordinates of x
 * @param yvec vector<double>:  coordinates of y
 * @param period not used, needed for consistency with Causality
 * @return bool : x-y timelike?
 */
bool Spacetime::Flat_causal(const vector<double>& xvec, const vector<double>& yvec,
                            vector<double>period, double mass)
{
    double dt = (xvec[0]-yvec[0]);
    double dspace2 = 0;
    for(int i=1; i<xvec.size(); i++){
       dspace2 += (xvec[i]-yvec[i])*(xvec[i]-yvec[i]);
    }
    if ((dt*dt) - (dspace2)  >0){
        return true;}
    else{
        return false;}
}

/**
 * @brief Function of two events in any D returning true if x-y timelike?,
 * else false in a periodic flat spacetime.
 * 
 * @param xvec : vector<double>, coordinates of x
 * @param yvec : vector<double>, coordinates of y
 * @param period : vector<double>, ith entry is periodicty along ith SPATIAL
 * dimension (0->x, 1->y, 2-->z)
 * @return bool : x-y timelike?
 */
bool Spacetime::Flat_causal_periodic(const vector<double>& xvec, 
                                        const vector<double>& yvec,
                                        vector<double> period,
                                        double mass)
{
    //std::cout<<"For debuggin: calling Flat_causal_periodic\n";
    double t_delta = (yvec[0] - xvec[0]);
    double t_delta2 = t_delta*t_delta;
    double space_delta2 = 0;
    for (int i = 1; i<yvec.size(); i++)
    {
        double space_delta_i = (yvec[i]-xvec[i]);
        double wrapped_space_delta_i = (space_delta_i>0)?
                                        period[i]-space_delta_i:
                                        period[i]+space_delta_i; 
        space_delta2 += (space_delta_i<wrapped_space_delta_i)?
                        space_delta_i*space_delta_i :
                        wrapped_space_delta_i*wrapped_space_delta_i;
    }
    return t_delta2 >= space_delta2 ;
}



/**
 * @brief Function of two events in any D returning {x-y timelike?, x<=y, x>y}.
 * 
 * @param xvec vector<double>:  coordinates of x
 * @param yvec vector<double>:  coordinates of y
 * @param period not used, needed for consistency with Causality
 * @return vector<bool> : {x-y timelike, x<=y, x>y}
 */
vector<bool> Spacetime::Flat_general_causal(const vector<double>& xvec, 
                                            const vector<double>& yvec,
                                            vector<double>period, double mass)
{
    double dt2 = (xvec[0]-yvec[0])*(xvec[0]-yvec[0]);
    double dspace2 = 0;
    for(int i=1; i<4; i++){
        dspace2 += (xvec[i]-yvec[i])*(xvec[i]-yvec[i]);
    }
    if (dt2-dspace2>0)
    {
        if (xvec[0]<yvec[0]) return {true,true, false};
        else return {true, false, true};
    }
    else
    {
        return {false, false, false};
    }
}


/**
 * @brief Function of two events in any D returning {x-y timelike?, x<=y, x>y}.
 * 
 * @param xvec : vector<double>, coordinates of x
 * @param yvec : vector<double>, coordinates of y
 * @param period : vector<double>, ith entry is periodicty along ith SPATIAL
 * dimension (0->x, 1->y, 2-->z)
 * @return vector<bool> : {x-y timelike, x<=y, x>y}
 */
vector<bool> Spacetime::Flat_general_causal_periodic(const vector<double>& xvec, 
                                                    const vector<double>& yvec,
                                                    vector<double> period,
                                                    double mass)
{
    //std::cout<<"For debuggin: calling Flat_causal_periodic\n";
    double t_delta = (yvec[0] - xvec[0]);
    double t_delta2 = t_delta*t_delta;
    double space_delta2 = 0;
    for (int i = 1; i<yvec.size(); i++)
    {
        double space_delta_i = (yvec[i]-xvec[i]);
        double wrapped_space_delta_i = (space_delta_i>0)?
                                        period[i]-space_delta_i:
                                        period[i]+space_delta_i; 
        space_delta2 += (space_delta_i<wrapped_space_delta_i)?
                        space_delta_i*space_delta_i :
                        wrapped_space_delta_i*wrapped_space_delta_i;
    }
    bool isCausal = t_delta2 >= space_delta2;
    return {isCausal,
            (t_delta >= 0.0) && isCausal,
            (t_delta < 0.0) && isCausal
           };
}