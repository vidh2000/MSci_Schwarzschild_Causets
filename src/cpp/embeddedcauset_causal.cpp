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
#include <string.h>
#include <vector>
#include <chrono>
#include <unordered_set>
#include <iterator>
#include <omp.h>
#include <stdint.h>
#include <utility>

#include "functions.h"
#include "vecfunctions.h"
#include "causet.h"
#include "embeddedcauset.h"
#include "shapes.h"
#include "spacetimes.h"

using std::vector;
using std::set;
using std::unordered_set;


/**
 * @brief Causal relation according to the spacetime.
 * 
 * @param xvec: vector<double>, coordinates of x
 * @param yvec: vector<double>, coordinates of y
 * 
 * @return : vector<bool> {x-y timelike, x<=y, x>y}.
 * @exception: returned if size of xvec and yvec diffrent than dimension of
 * spacetime.
 */
bool EmbeddedCauset::causality(vector<double> xvec, 
                                vector<double> yvec)
{
    auto xycausality = this->_spacetime.Causality();
    return xycausality(xvec, yvec, _spacetime._period, _spacetime._mass);
};


/**
 * @brief Causal relation according to the spacetime.
 * 
 * @param xvec: vector<double>, coordinates of x
 * @param yvec: vector<double>, coordinates of y
 * 
 * @return : vector<bool> {x-y timelike, x<=y, x>y}.
 * @exception: returned if size of xvec and yvec diffrent than dimension of
 * spacetime.
 */
std::vector<bool> EmbeddedCauset::general_causality(vector<double> xvec, 
                                                    vector<double> yvec)
{
    auto xycausality = this->_spacetime.General_Causality();
    return xycausality(xvec, yvec, _spacetime._period, _spacetime._mass);
};



/**
 * @brief Causal relation according to the spacetime.
 * 
 * @param xvec: vector<double>, coordinates of x
 * @param yvec: vector<double>, coordinates of y
 * 
 * @return Bool: true if two events are timelike, else false.
 * @exception: returned if size of xvec and yvec difefrent than dimension of
 * spacetime.
 */
bool EmbeddedCauset::areTimelike4D(vector<double> &xvec, vector<double> &yvec)
{
    double dt = (xvec[0]-yvec[0]);
    double dspacex = xvec[1]-yvec[1];
    double dspacey = xvec[2]-yvec[2];
    double dspacez = xvec[3]-yvec[3];
    // for(int i=1; i<dim; i++){
    //    dspace2 += (xvec[i]-yvec[i])*(xvec[i]-yvec[i]);
    // }
    if ((dt*dt)-(dspacex*dspacex + dspacey*dspacey + dspacez*dspacez)>0)
    {
        return true;}
    else{
        return false;}
};


/**
 * @brief Causal relation according to the spacetime.
 * 
 * @param xvec: vector<double>, coordinates of x
 * @param yvec: vector<double>, coordinates of y
 * 
 * @return Bool: true if event x preceeds y.
 * @exception: returned if size of xvec and yvec difefrent than dimension of
 * spacetime.
 */
bool EmbeddedCauset::AprecB(vector<double> xvec, vector<double> yvec)
{
    auto atimelikeb = _spacetime.General_Causality();
    return atimelikeb(xvec, yvec, _spacetime._period, _spacetime._mass)[1];
};