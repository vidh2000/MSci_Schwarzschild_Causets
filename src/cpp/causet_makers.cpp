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


void Causet::make_sets_fromC()
{
    _futures.resize(_size);
    _pasts.resize(_size);
    // Loop through coordinates t_min -> t_max.
    // j>i automatically imposed as C_ij <-> i precedes j. 
    for (int j = 1; j<_size; j++)
    {
        for(int i=j-1; i>-1; i--)
        {
            if (_CMatrix[i][j]!=0) 
            {
                // Add i and its past to the past of j
                _pasts[j].insert(i);
                _futures[i].insert(j);
            }
        }
    }
}

void Causet::make_cmatrix(){}
void Causet::make_pasts(){}
void Causet::make_futures(){}
void Causet::make_past_links(){}

void Causet::make_future_links_fromC()
{
    if (_CMatrix.size()==0)
    {
        std::cout << "To create future link matrix, CMatrix must exist";
        throw std::invalid_argument("No CMatrix");}

    _future_links.resize(_size);
    
    #pragma omp parallel for
    for (int i=0; i<_size; i++)
    {
        for (int j=i+1; j<_size; j++)
        {
            if (_CMatrix[i][j] == 0) {
                continue;
            }
            else
            {
                bool has_broken = false;
                for (int k=i+1; k<j;k++)
                {
                    if (_CMatrix[i][k]*_CMatrix[k][j]!=0){
                        has_broken = true;
                        break;}
                }
                if (!has_broken){
                    _future_links[i].insert(j);}
            }
        }
    }
}