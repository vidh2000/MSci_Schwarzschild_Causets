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
#include <random>
#include <stdint.h>
#include "functions.h"
#include "vecfunctions.h"
#include "causet.h"

using std::vector;
using std::set;
using std::unordered_set;

/**
 * @brief Yields a copy of the interval's cmatrix obtained by removing
 *          labels of elements not belonging to said interval.
 * 
 * @param ordered_interval Interval (vector) of elements that remain
 * @param make_matrix, make_sets, make_links are booleans saying
 *          what type of causet was created (what exists -> to know what
 *                                                          to update)
 */
std::vector<vector<int>> Causet::getIntervalCmatrix(
                                    vector<int> ordered_interval)
{
    if (_CMatrix.size())
    {
        return get_reducedMatrix(_CMatrix, ordered_interval);
    } 
    else 
    {
        print("Require existing CMatrix! It doesn't exist.");
        throw std::runtime_error("");
    }
}


/**
 * @brief Given an element x of the causet, compute N_k(x), i.e. the number of
 * y in causet such that |I(y,x)|= k + 1. Return an array of N_k(x) for 
 * k in [kmin, kmax].
 * 
 * @param x int : label of element with respect to which find size of layers.
 * @param kmax int : k of the maximum k-past-layer to count.
 * @param kmin int : k of the minimum k-past-layer to count. Default is 1, which
 * corresponds to the count of past links.
 * 
 * @return std::vector<double> Nk_BD : vector of Nk_s(x) for k in [kmin, kmax].
 */
std::vector<double> Causet::Nk_BD (int x, int kmax, int kmin)
{
    // Define vector where to store values
    std::vector<double> Nk_s (kmax-kmin+1);

    // Loop over elements y potentially in past of x
    for (int y = x-1; y>-1; y--)
    {
        // if i prec x, get size of interval
        // if size in [kmin, kmax], update count
        if (_CMatrix[y][x] != 0)
        {
            int n_yx = IntervalCard(y, x);
            int k_yx = n_yx - 1;
            if (kmin <= k_yx && k_yx <= kmax)
            {
                Nk_s[k_yx] += 1;
            }
        } 
    }

    return Nk_s;
}


/**
 * @brief Compute cardinality of causality interval between a and b.
 * 
 * @param a int : label of event a.
 * @param b int : label of event b
 * @param includeBoundary bool : innclude a and b in count?
 * @return int : Cardinality of Interval. 0 if they are not related.
 */
int Causet::IntervalCard(int a, int b, bool includeBoundary)
{
    if (a==b)
        {return 1;}
    // IF DEFINED USE CMATRIX
    if (_CMatrix.size())
    {
        if ((a<b && _CMatrix[a][b] == 0) || (b<a && _CMatrix[b][a] == 0)){
            return 0;
        }

        int Nintersections = 2 * includeBoundary;
        if (a<b)
        {
            for (int i = a+1; i<b; i++)
                {Nintersections += _CMatrix[a][i] && _CMatrix[i][b];}
        }
        else
        {
            for (int i = b+1; i<a; i++)
                {Nintersections += _CMatrix[i][a] && _CMatrix[b][i];}
        }
        return Nintersections;
    }
    // ELSE USE SETS
    else //if ( _pasts.size() && _futures.size()) and no cmatrix exists
    {

        if (a<b)
        {
            if (_futures[a].find(b)!=_futures[a].end()){
                return 0;
            }
            int Nintersections = 2 * includeBoundary;
            if (_futures[a].size()<_pasts[b].size()) //loop over shortest
            {
                for (int e_ai : _futures[a])
                    {Nintersections += _pasts[b].find(e_ai) !=_pasts[b].end();}
            }
            else
            {
                for (int e_bi : _pasts[b])
                    {Nintersections+= _futures[a].find(e_bi)!=_futures[a].end();}
            }
            return Nintersections;
        }

        else /*b<a*/
        {
            if ( _futures[b].find(a)!=_futures[b].end()){
                return 0;
            }
            int Nintersections = 2 * includeBoundary;
            if (_pasts[a].size()<_futures[b].size()) //loop over shortest
            {
                for (int e_ai : _pasts[a])
                    {Nintersections+= _futures[b].find(e_ai)!=_futures[b].end();}
            }
            else
            {
                for (int e_bi : _futures[b])
                    {Nintersections += _pasts[a].find(e_bi) !=_pasts[a].end();}
            }
            return Nintersections;
        }
    }  
}