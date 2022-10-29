#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <stdio.h>
#include <stack>
#include <string>
#include <vector>

#include "shapes.h"
#include "spacetimes.h"
#include "functions.h"
#include "MyVecFunctions.h"

/*============================================================================
* ============================================================================
* GENERAL SPACETIME METHODS
* ============================================================================
============================================================================*/
double Spacetime::Parameter (char key)
/**
 * @brief Return parameter "key" for the shape of the spacetime
 */
    {return _params[key];}

CoordinateShape Spacetime::DefaultShape()
/**
 * @brief Returns default coordinate shape of embedding region in spcetime.  
 */
    {return CoordinateShape(Dim, 'cylinder');}


std::vector<double> Spacetime::_T_slice_sampling(double t, 
                                                std::vector<double>origin,
                                                int samplingsize = 128)
/**
 * @brief Internal function for the time sampling array for a cone from 
 * "origin" to time "t"
 */
    {return linspace(origin[0], t, samplingsize);}     


std::vector<bool> Spacetime::causal1d(std::vector<double> xvec, 
                                      std::vector<double> yvec,
                                      double causality_eps)
/**
 * @brief The most simple causal function: t_2>t_1
 */
{
    double t_delta = yvec[0] - xvec[0];
    return {t_delta >= 0.0, t_delta < 0.0};
}

typedef std::vector<bool> (*func)
(std::vector<double> xvec, std::vector<double> yvec);
func Spacetime::Causality()
/**
 * @brief Return most simple causal function: t_2>t_1, useful
 * in 1D cases. 
 */
    {return &Spacetime::causal1d;}




/*============================================================================
* ============================================================================
* FLAT SPACETIME
* Notes: does not support period
* ============================================================================
============================================================================*/
FlatSpacetime::FlatSpacetime(int dim = 2)
{
    if (Dim < 1)
    {
        throw std::invalid_argument("Dimension has to be at least 1.");
    }

    _dim = dim;
    _name = 'flat';
    _metricname = 'Minkowski';
}

CoordinateShape FlatSpacetime::DefaultShape()
    {return CoordinateShape(_dim, 'diamond');}


double FlatSpacetime::ds2(std::vector<double> xvec, std::vector<double> yvec)
/**
 * @brief Spacetime Interval Squared
 */
{
    
    double time_delta2  = (yvec[0] - xvec[0])*(yvec[0] - xvec[0]);
    double space_delta2 = 0;
    for (int i = 1; i<_dim; i++)
        {space_delta2 += (yvec[i]-xvec[i])*(yvec[i]-xvec[i]);}
    return time_delta2 - space_delta2;
}

double FlatSpacetime::ds(std::vector<double> xvec, std::vector<double> yvec)
/**
 * @brief Spacetime Interval
 */
{
    
    double time_delta2  = (yvec[0] - xvec[0])*(yvec[0] - xvec[0]);
    double space_delta2 = 0;
    for (int i = 1; i<_dim; i++)
        {space_delta2 += (yvec[i]-xvec[i])*(yvec[i]-xvec[i]);}
    return std::sqrt(time_delta2 - space_delta2);
}

     
typedef std::vector<bool> (*func)
(std::vector<double> xvec, std::vector<double> yvec, double causality_eps);
func FlatSpacetime::Causality()
/**
 * @brief Return {bool1, bool2} &callable (vector<double> xvec, yvec)
 *        which takes coordinates of 2 points and returns vector with booleans
 *        {x<=y, x>y}
 */
{
    if (_dim == 1)
        {return &Spacetime::causal1d;}
    else
        {return &FlatSpacetime::causal;}
}

std::vector<bool> FlatSpacetime::causal (std::vector<double> xvec, 
                                         std::vector<double> yvec,
                                         double causality_eps = 1e-12)
/**
 * @brief Function of two events in any D returning {x<=y, x>y} with error in 
 *        calculations causality_eps (1e-12 by default)
 */
{
    double t_delta = (yvec[0] - xvec[0]);
    double t_delta2 = t_delta*t_delta;
    double space_delta2 = 0;
    for (int i = 1; i<yvec.size(); i++)
        {space_delta2 += (yvec[i]-xvec[i])*(yvec[i]-xvec[i]);}
    bool isCausal = t_delta2 >= space_delta2 - causality_eps;
    return {(t_delta >= 0.0) && isCausal, (t_delta < 0.0) && isCausal};
}




/*============================================================================
* ============================================================================
* BLACK HOLES SPACETIME
* Notes: does not support anything actually at the moment
* ============================================================================
============================================================================*/
BlackHoleSpacetime::BlackHoleSpacetime(int dim = 2, double r_S = 0.5,
                                       char metric = 'Eddington-Finkelstein')
/**
 * @brief Initialises a Black Hole Spacetime.
 *
 * @param dim: dimension of spacetime. Default 2.
 * @param r_S: Schwarzschild radius. Default 0.5
 * @param metric: Specify metric: either 'Eddington-Finkelstein' or 'EF'
 *                (default), or 'Schwarzschild' or 'S'.
 */
{
    if (Dim < 2 || Dim > 4)
    {
        throw std::invalid_argument("Dimension has to be 2, 3 or 4.");
    }
    _dim = dim;
    _name = 'black hole';

    if (metric == 'Eddington-Finkelstein' || metric == 'EF')
        {_metricname = 'Eddington-Finkelstein';}
    else if (metric == 'Schwarzschild' || metric == 'S')
        {_metricname = 'Schwarzschild';}
}
