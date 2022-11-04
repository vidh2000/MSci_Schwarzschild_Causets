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
#include <vector>
#include <chrono>
#include <unordered_set>

//using namespace std::chrono;


#include "functions.h"
#include "MyVecFunctions.h"
#include "causet.h"
#include "embeddedcauset.h"
#include "shapes.h"
#include "spacetimes.h"

using std::vector;
using std::set;
using std::unordered_set;


//=============================================================================
//CONSTRUCTORS  //=============================================================
//=============================================================================



//=============================================================================
//MAKE ATTRS  //===============================================================
//=============================================================================
/**
 * @brief Creates _pasts: the set of past events for each event
 * 
 * @param method: const char*, possible choices are
 * - "coordinates": create from coordinates causality
 * - "Cmatrix": create from already existing _CMatrix
 * - "futures": create from already existing futures
 */
void EmbeddedCauset::make_all_pasts(const char* method = "coordinates")
{   
    vector<std::unordered_set<int>> pasts;
    pasts.resize(_size);
    past_links.resize(size);
    //std::cout << "Finished resizing sets.." << std::endl;
    // Loop through coordinates t_min -> t_max
    for (int i=1; i<size; i++)
    {
        //std::cout << "Event #"<< i+1 << std::endl;
        for(int j=i-1; j>-1; j--)
        {
            // Check if j^th element is in pasts[i]
            //if (set_contains(j, pasts[i])){
            //    //std::cout << "For i=" << i << "j=" << j << "is contained"
            //    //                    << std::endl;
            //    continue;
            //}
            //else{
            // Check if j<i (are causally connected)
            if (areTimelike(coords[i],coords[j])){
                // If yes, add e_j and its past to the past of e_i
                pasts[i].insert(j);
                pasts[i].insert(pasts[j].begin(), pasts[j].end());
            }
            //}
        }
    }
    std::cout <<"Finished sprinkling..." << std::endl;
}


void EmbeddedCauset::make_cmatrix()
{

    std::cout << "Creating causet N=" << size << " via cmatrix" << std::endl;

    // Creating a 
    vector<vector<int>> cmatrix;
    cmatrix.resize(size);
    for(int i=0; i<size; i++) {
        cmatrix[i].resize(size);
    }
    for (int i=1; i<size; i++)
    {
        //std::cout << "Event #"<< i+1 << std::endl;
        for(int j=0; j<i; j++){
            if (areTimelike(coords[i],coords[j]))
            {
                cmatrix[i][j] = 1;
            }
        }
    }
}


//=============================================================================
//GETTERS     //===============================================================
//=============================================================================


//=============================================================================
//RELATIONS   //===============================================================
//=============================================================================
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
bool EmbeddedCauset::areTimelike(vector<double> xvec, vector<double> yvec)
{
    auto atimelikeb = _spacetime.Causality();
    return atimelikeb(xvec, yvec)[0];
    // int dim = xvec.size();
    // double time_delta2  = (yvec[0] - xvec[0])*(yvec[0] - xvec[0]);
    // double space_delta2 = 0;
    // for (int i = 1; i<dim; i++){
    //     space_delta2 += (yvec[i]-xvec[i])*(yvec[i]-xvec[i]);
    // }
    // if (time_delta2 - space_delta2 >0){
    //     return true;
    // }
    // else{
    //     return false;
    // }
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
    auto atimelikeb = _spacetime.Causality();
    return atimelikeb(xvec, yvec)[1];
};


//=============================================================================
//MODIFY      //===============================================================
//=============================================================================

