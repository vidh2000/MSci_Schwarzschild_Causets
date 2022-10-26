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
#include "spacetimes.h"

/*============================================================================
* GENERAL SPACETIME METHODS
============================================================================*/
float Spacetime::Parameter (char key)
/**
 * @brief Return parameter "key" for the shape of the spacetime
 */
{
    return _params[key];
}

CoordinateShape Spacetime::DefaultShape()
/**
 * @brief Returns default coordinate shape of embedding region in spcetime.  
 */
{
    return CoordinateShape(Dim, 'cylinder');
}

std::vector<float> Spacetime::_T_slice_sampling(float t, 
                                                std::vector<float>origin,
                                                int samplingsize = 128)
/**
 * @brief Internal function for the time sampling array for a cone from 
 * "origin" to time "t"
 */
{
    ;
}     
