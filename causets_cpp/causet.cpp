/// \authors Vid Homsak, Stefano Veroni
/// \date 29/09/2022

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

// For Vid (can't locate headers for some reason)
// Path: D:\Documents\Sola\Imperial College London\Year 4\MSci project\Project\causets_code\causets_cpp\"header".h...

#include "functions.h"
#include "MyVecFunctions.h"
#include "causet.h"

// #include <D:\Documents\Sola\Imperial College London\Year 4\MSci project\Project\causets_code\causets_cpp\functions.h>
// #include <D:\Documents\Sola\Imperial College London\Year 4\MSci project\Project\causets_code\causets_cpp\causet.h>

using std::vector;
using std::set;


//============================================================================
//////////////////////////////////////////////////////////////////////////////
// CONSTRUCTORS   ////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//============================================================================

/**
 * @brief Default constructor.
 * 
 * @param future : vector<set<int>. ith entry is set of events at future of i.
 * @param past :  vector<set<int>. ith entry is set of events at past of i.
 */
Causet::Causet(vector<set<int>> future = {}, vector<set<int>> past = {})
    {
        _past   = past;
        _future = future;
    }


/**
 * @brief Construct a new Causet object from Causal (or Link) matrix.\n
 *        The Causal matrix is C_ij = 1 iff j in i's Future.\n
 *        The Link matrix is L_ij = 1 iff j in i's Future and ij is a link.
 * 
 * @param M : vector<vector<int>>. The matrix.
 * @param Lmatrix : bool. Default false.
 *                  If true (Not implemented) construct from Link matrix.
 * @param reverse_causality : bool. Default false.
 *                  If true generate from trasnpose of Causal (or link) matrix.
 * 
 * @exception invalid_argument thrown is M is not a square matrix.
 * @exception invalid_argument thrown is M is not anti-symmetric.
 */
Causet::Causet(vector<vector<int>> M, bool Lmatrix = false, 
                bool reverse_causality = false)
{
    int rowcount = M.size();
    int colcount = M[0].size();
    if (colcount != rowcount)
    {
        throw std::invalid_argument
        ("The given matrix is not a square matrix!");
    }
    
    vector<set<int>> _future (colcount);
    vector<set<int>> _past   (colcount);

    for (int i = 0; i<colcount; i++)
    {
        for (int j = i+1; j<colcount; j++)
        {
            if (M[i][j] && M[j][i])
            {
                throw std::invalid_argument
                ("The given matrix is not antisymmetric!");
            }
            else if (M[i][j])
            {
                _future[i].insert(j);
                _past[j].insert(i);
            }
            else if (M[j][i])
            {
                _future[j].insert(i);
                _past[i].insert(j);
            }
        }
    }
}

/**
 * @brief Construct a new Causet object from Causal (or Link) matrix.\n
 *        The Causal matrix is C_ij = 1 iff j in i's Future.\n
 *        The Link matrix is L_ij = 1 iff j in i's Future and ij is a link.
 * 
 * @param M : vector<vector<bool>>. The matrix. 
 * @param Lmatrix : bool. Default false.
 *                  If true (Not implemented) construct from Link matrix.
 * @param reverse_causality : bool. Default false.
 *                  If true generate from trasnpose of Causal (or link) matrix.
 * 
 * @exception invalid_argument thrown is M is not a square matrix.
 * @exception invalid_argument thrown is M is not anti-symmetric.
 */
Causet::Causet(vector<vector<bool>> M, bool Lmatrix = false, 
                bool reverse_causality = false)
{
    int rowcount = M.size();
    int colcount = M[0].size();
    if (colcount != rowcount)
    {
        throw std::invalid_argument
        ("The given matrix is not a square matrix!");
    }
    
    vector<set<int>> _future (colcount);
    vector<set<int>> _past   (colcount);

    for (int i = 0; i<colcount; i++)
    {
        for (int j = i+1; j<colcount; j++)
        {
            if (M[i][j] && M[j][i])
            {
                throw std::invalid_argument
                ("The given matrix is not antisymmetric!");
            }
            else if (M[i][j])
            {
                _future[i].insert(j);
                _past[j].insert(i);
            }
            else if (M[j][i])
            {
                _future[j].insert(i);
                _past[i].insert(j);
            }
        }
    }
}

//============================================================================
//////////////////////////////////////////////////////////////////////////////
// MODIFICATIONS  ////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//============================================================================








