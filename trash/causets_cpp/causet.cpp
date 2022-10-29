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
#include "causetevent.h"
#include "causet.h"

// #include <D:\Documents\Sola\Imperial College London\Year 4\MSci project\Project\causets_code\causets_cpp\functions.h>
// #include <D:\Documents\Sola\Imperial College London\Year 4\MSci project\Project\causets_code\causets_cpp\causet.h>


//////////////////////////////////////////////////////////////////////////////
// CONSTRUCTORS   ////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

/**
 * @brief Default constructor from set of causet events
*/
Causet::Causet(set<CausetEvent> set = {})
    {_events = set;}


/**
 * @brief Construct a new Causet object from Causal (or Link) matrix.\n
 *        The Causal matrix is C_ij = 1 iff j in i's Future.\n
 *        The Link matrix is L_ij = 1 iff j in i's Future and ij is a link.
 * 
 * @param Lmatrix : bool. Default false.
 *                  If true (Not implemented) construct from Link matrix.
 * @param reverse_causality : bool. Default false.
 *                  If true generate from trasnpose of Causal (or link) matrix.
 * 
 * @exception invalid_argument thrown is M is not a square matrix.
 * @exception invalid_argument thrown is M is not anti-symmetric.
 */
Causet::Causet(vector<vector<double>> M, bool Lmatrix = false, 
                bool reverse_causality = false)
{
    int rowcount = M.size();
    int colcount = M[0].size();
    if (colcount != rowcount)
        {throw invalid_argument("The given matrix is not a square matrix!");}
    
    set<CausetEvent> events;
    for (int i = 0; i<colcount; i++)
    {
        CausetEvent e = CausetEvent(i);
        events.insert(e);
    }
    for (int i = 0; i<colcount; i++)
    {
        CausetEvent e = events[i];
        for (int j = i+1; j<colcount; j++)
        {
            if (M[i][j] == 1 && M[j][i] == 1)
                {throw invalid_argument("The given matrix is not \
                antisymmetric!");}
            else if (M[i][j])
                {e._addToFuture(events[j]);}
        }
    }
}








