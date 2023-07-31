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
 * @brief Construct a new Causet:: Causet object
 * 
 * @param Cmatrix 
 */
Causet::Causet(vector<vector<int>> Cmatrix)
{
    _CMatrix = Cmatrix;
    _size = Cmatrix.size();
}

/**
 * @brief Construct Causet Object.
 * 
 * @param Cmatrix 
 */
template <typename num>
Causet::Causet(vector<vector<num>> Cmatrix)
{
    _CMatrix = (int)Cmatrix;
    _size = Cmatrix.size();
}


std::vector<std::vector<int>> Causet::get_CMatrix(){
    return _CMatrix;
}

std::vector<std::vector<int>> Causet::CMatrix(std::vector<int> labels)
{
    if (!labels.size())
        {return _CMatrix;}
    else
        {return {};}
}


int Causet::size() 
    {return _size;}
bool Causet::is_CMatrix_special()
    {return _special_matrix;}
bool Causet::is_Cij_special()
    {return _special_matrix;}


void Causet::saveC(const char* path_file_ext)
{
    std::ofstream out;
    out.open(path_file_ext);
    for (auto row : _CMatrix) 
    {
        for (auto col : row)
            {out << col <<',';}
        out<<std::endl;
    }
    out.close();
    return;
}