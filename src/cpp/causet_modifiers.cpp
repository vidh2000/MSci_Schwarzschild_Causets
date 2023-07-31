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
 * @brief Coarse grain Causet of "card" events.
 * 
 * @param card : int
 * @param make_matrix : if true and _CMatrix non-empty, update _CMatrix 
 * @param make_sets : if true, update _pasts and/or _futures, if they are 
 * defined
 * @param make_links : if true, update _past and/or _future links, if they are 
 * defined
 */
void Causet::coarsegrain(int card, bool make_matrix, 
                         bool make_sets, bool make_links, int seed)
{
    vector<int> labels = distinct_randint(card, _size, seed);
    this->discard(labels, make_matrix, make_sets, make_links); 
}
void Causet::cgrain(int card, bool make_matrix, 
                    bool make_sets, bool make_links, int seed)
{
    vector<int> labels = distinct_randint(card, _size, seed);
    this->discard(labels, make_matrix, make_sets, make_links); 
}

/**
 * @brief Coarse grain Causet of "fract*_size" events.
 * 
 * @param fract : double, fraction of elements to remove
 * @param make_matrix bool: if true and _CMatrix non-empty, update _CMatrix 
 * @param make_sets bool: if true, update _pasts and/or _futures, if they are 
 * defined
 * @param make_links bool: if true, update _past and/or _future links, if they 
 * are defined
 */
void Causet::coarsegrain(double fract, bool make_matrix, 
                         bool make_sets, bool make_links, int seed)
{
    int card = fract * _size;
    vector<int> labels = distinct_randint(card, _size, seed);
    this->discard(labels, make_matrix, make_sets, make_links); 
}
void Causet::cgrain(double fract, bool make_matrix, 
                    bool make_sets, bool make_links, int seed)
{
    int card = fract * _size;
    vector<int> labels = distinct_randint(card, _size, seed);
    this->discard(labels, make_matrix, make_sets, make_links); 
}

/**
 * @brief Discard event number "label" (from 0 to N-1)
 * 
 * @param label : int
 * @param make_matrix : if true, update _CMatrix, if already defned
 * @param make_sets : if true, update _pasts and/or _futures, if they are 
 * defined
 * @param make_links : if true, update _past and/or _future links, if they are 
 * defined
 */
void Causet::discard(int label, bool make_matrix, 
                     bool make_sets, bool make_links)
{
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
    _dim = 0;
}


/**
 * @brief Discard events numbered by i in "label" (from 0 to N-1)
 * 
 * @param label : vector<int>
 * @param make_matrix : if true, update _CMatrix, if already defned
 * @param make_sets : if true, update _pasts and/or _futures, if they are 
 * defined
 * @param make_links : if true, update _past and/or _future links, if they are 
 * defined
 */
void Causet::discard(vector<int> labels, bool make_matrix, 
                    bool make_sets, bool make_links)
{
    if (make_matrix)
    {
        if (_CMatrix.size())
        {
            remove_indices(_CMatrix, labels);
            for (vector<int> row : _CMatrix)
                {remove_indices(row, labels);}
        } 
    }
    if (make_sets)
    {
        if (_pasts.size())
        {
            remove_indices(_pasts, labels);
            for (unordered_set<int> past_i : _pasts)
                {discard_from_set(past_i, labels);}
        } 
        if (_futures.size())
        {
            remove_indices(_futures, labels);
            for (unordered_set<int> fut_i : _futures)
                {discard_from_set(fut_i, labels);}
        }   
    }
    if (make_links)
    {
        if (_past_links.size())
        {
            remove_indices(_past_links, labels);
            for (unordered_set<int> plinks_i : _past_links)
                {discard_from_set(plinks_i, labels);}
        } 
        if (_future_links.size())
        {
            remove_indices(_future_links, labels);
            for (unordered_set<int> flinks_i : _future_links)
                {discard_from_set(flinks_i, labels);}
        }   
    }
    _size--;
    _dim = 0;
} 