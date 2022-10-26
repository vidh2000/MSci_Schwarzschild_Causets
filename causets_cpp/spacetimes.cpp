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
 * 
 */
{

};
virtual CoordinateShape DefaultShape();

typedef bool (*func)(std::vector<float> xvec, std::vector<float> yvec);
virtual func  Causality();   

virtual std::vector<float>T_slice_sampling(float t, 
                                            std::vector<float>origin,
                                            int samplingsize = -1);     
};