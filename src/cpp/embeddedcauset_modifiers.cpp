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
 * @brief Sort coordinates of Causet.
 * 
 * @param dim : int. Dimension to use for the sorting.
 * @param reverse : bool. If true (non-default) reverse order.
 * @return vector<vector<double>> 
 */
void EmbeddedCauset::sort_coords(int dim,// = 0,
                                 bool reverse)// = false)
{
    int rev_factor = (reverse)? -1 : 1;
    auto sort_lambda = [dim, rev_factor] (vector<double> v1, vector<double> v2) 
                      {return rev_factor * v1[dim]<v2[dim];};
    std::sort(_coords.begin(), _coords.end(), sort_lambda);
    //std::cout << "Coords after sorting:" << std::endl;
    //print_vector(_coords);
    //std::cout << std::endl;
}


void EmbeddedCauset::add(vector<double> xvec)
{
    // Update CMatrix (if defined)
    // Remove from sets and scale all following one down (if defined)
    // (maybe redefining sets is faster)
    // Increase size by one
    // Turn _dim to 0 as new causet
}


void EmbeddedCauset::discard(int label, bool make_matrix, // = true, 
                             bool make_sets, // = false,
                             bool make_links) // = true)
{
    print("Don't think this one works tbf, needs to be tested if needed");
    _coords.erase(_coords.begin() + label);

    if (make_matrix)
    {
        if (_CMatrix.size())
        {
            _CMatrix.erase(_CMatrix.begin() + label);
            for (vector<int> row : _CMatrix)
                {row.erase(row.begin() + label);}
        } 
    }
    if (make_sets)
    {
        if (_pasts.size())
        {
            _pasts.erase(_pasts.begin()+label);
            for (unordered_set<int> past_i : _pasts)
                {discard_from_set(past_i, label);}
        } 
        if (_futures.size())
        {
            _futures.erase(_futures.begin()+label);
            for (unordered_set<int> fut_i : _futures)
                {discard_from_set(fut_i, label);}
        }   
    }
    if (make_links)
    {
        if (_past_links.size())
        {
            _past_links.erase(_past_links.begin()+label);
            for (unordered_set<int> plinks_i : _past_links)
                {discard_from_set(plinks_i, label);}
        } 
        if (_future_links.size())
        {
            _future_links.erase(_future_links.begin()+label);
            for (unordered_set<int> flinks_i : _future_links)
                {discard_from_set(flinks_i, label);}
        }   
    }
    _size--;
    //_dim = 0;
}


/**
 * @brief Cuts the cmatrix and relabels the futures(links) and pasts(links)
 *        to reduce the whole causet to a smaller interval/removes certain
 *        labels
 * 
 * @param labels Discard labels (get rid of this? not used now)
 * @param ordered_interval Interval (vector) of elements that remain
 * @param make_matrix, make_sets, make_links are booleans saying
 *          what type of causet was created (what exists -> to know what
 *                                                          to update)
 */
void EmbeddedCauset::discard(vector<int> labels, //get rid of this?
                             vector<int> ordered_interval,
                             bool make_matrix, // = true, 
                             bool make_sets, // = false,
                             bool make_links) // = true)
{
    remove_indices(_coords, labels);


    if (make_matrix)
    {
        if (_CMatrix.size())
        {
            remove_indices_fromCmatrix(_CMatrix, ordered_interval);
            //for (vector<int> row : _CMatrix)
            //    {remove_indices(row, labels);}
        } 
    }
    if (make_sets)
    {
        if (_pasts.size())
        {
            replace_indices(_pasts, ordered_interval);
            //for (unordered_set<int> past_i : _pasts)
            //    {discard_from_set(past_i, labels);}
        } 
        if (_futures.size())
        {
            replace_indices(_futures, ordered_interval);
            //for (unordered_set<int> fut_i : _futures)
            //    {discard_from_set(fut_i, labels);}
        }   
    }
    if (make_links)
    {
        print("NOTE: Discarding links is not yet checked (but updated)!");
        if (_past_links.size())
        {
            replace_indices(_past_links, ordered_interval);
            //for (unordered_set<int> plinks_i : _past_links)
            //    {discard_from_set(plinks_i, labels);}
        } 
        if (_future_links.size())
        {
            replace_indices(_future_links, ordered_interval);
            //for (unordered_set<int> flinks_i : _future_links)
            //    {discard_from_set(flinks_i, labels);}
        }   
    }

    _size -= labels.size();
    //_dim = 0;
}