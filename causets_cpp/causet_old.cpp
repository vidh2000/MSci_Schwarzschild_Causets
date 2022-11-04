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
 * @param pasts :  vector<set<int>. ith entry is set of events at past of i.
 */
Causet::Causet(vector<set<int>> pasts = {})
{
    _pasts   = pasts;
}

/**
 * @brief Constructor of both past and future.
 * 
 * @param pasts :  vector<set<int>. ith entry is set of events at past of i.
 * @param other : vector<set<int>. ith entry is set of events at future of i.
 * @param isfuture : bool. If true other is _future, otherwise _past_links.
 */
Causet::Causet(vector<set<int>> pasts, vector<set<int>> other,
                bool isfuture = true)
{
    _pasts   = pasts;
    if (isfuture)
        {_futures = other;}
    else
        {_past_links = other;}   
}

/**
 * @brief Constructor of both past and future and past_links.
 * 
 * @param pasts :  vector<set<int>. ith entry is set of events at past of i.
 * @param future : vector<set<int>. ith entry is set of events at future of i.
 * @param past_links: vector<set<int>>. ith entry is set of events at past of i
 * such that ij is a link.
 */
Causet::Causet(vector<set<int>> pasts, vector<set<int>> futures,
                vector<set<int>> past_links)
{
    _pasts   = pasts;
    _futures = futures;
}

/**
 * @brief Construct a new Causet object from Causal (or Link) matrix.\n
 *        The Causal matrix is C_ij = 1 iff j in i's Future.\n
 *        The Link matrix is L_ij = 1 iff j in i's Future and ij is a link.
 * 
 * @param M : vector<vector<int>>. The matrix.
 * @param Lmatrix : bool. Default false.
 *                  If true (Not implemented) construct from Link matrix.
 * @param reverse: bool. Default false.
 *                  If true generate from trasnpose of Causal (or link) matrix.
 * 
 * @exception invalid_argument thrown is M is not a square matrix.
 * @exception invalid_argument thrown is M is not anti-symmetric.
 */
Causet::Causet(vector<vector<int>> M, bool Lmatrix = false, 
                bool reverse = false)
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
            else if (M[i][j] || (M[j][i] && reverse))
            {
                _futures[i].insert(j);
                _pasts[j].insert(i);
            }
            else if (M[j][i] || (M[i][j] && reverse))
            {
                _futures[j].insert(i);
                _pasts[i].insert(j);
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
 * @param reverse: bool. Default false.
 *                  If true generate from trasnpose of Causal (or link) matrix.
 * 
 * @exception invalid_argument thrown is M is not a square matrix.
 * @exception invalid_argument thrown is M is not anti-symmetric.
 */
Causet::Causet(vector<vector<bool>> M, bool Lmatrix = false, 
                bool reverse = false)
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
            else if (M[i][j] || (M[j][i] && reverse))
            {
                _futures[i].insert(j);
                _pasts[j].insert(i);
            }
            else if (M[j][i] || (M[i][j] && reverse))
            {
                _futures[j].insert(i);
                _pasts[i].insert(j);
            }
        }
    }
}


//============================================================================
//////////////////////////////////////////////////////////////////////////////
// SETTERS   /////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//============================================================================

/**
 * @brief Update _futures attribute based on _pasts.
 * 
 */
void Causet::update_futures()
{
    vector<set<int>> futures (_pasts.size());
    for (int i = 0; i<_pasts.size(); i++)
    {
        set<int> p_i =_pasts[i];
        for (int j : p_i)
        {
            futures[j].insert(i);
        }
    }
    _futures = futures;
    return;
}

/**
 * @brief Update _past_links attribute based on _pasts and _futures. Compute
 * 
 * @param update_futures : bool. If True, update _futures first. _futures will
 *                         be updated anyways if empty.
 */
void Causet::update_links(bool update_futures = false)
{
    if (_futures.empty() || update_futures)
        {this->update_futures();}
    for (int i = 0; i<_pasts.size(); i++)
    {
        set<int> p_i =_pasts[i];
        for (int j : p_i)
        {
            set<int> inter_ij = {};
            std::set_intersection(p_i.begin(), p_i.end(),
                                _futures[j].begin(), _futures[j].end(),
                                std::inserter(inter_ij, inter_ij.end()));
            if (inter_ij.empty())
                {_past_links[i].insert(j);}
        }
    }
}


//============================================================================
//////////////////////////////////////////////////////////////////////////////
// MODIFICATIONS  ////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//============================================================================

/**
 * @brief Merge 2 causets with pastCauset in the past of futCauset. Note,
 *        futCauset's labels are shifted i --> i + pastCauset.size().
 * 
 * @param pastCauset : Causet. Goes in the past.
 * @param futCauset : Causet. Goes in the future.
 * @param disjoint : bool. Default false. If false connect the two causets, 
 *                  otherwise leave them separated.
 * @param update_futures : bool. Default false. If true update future sets of
 *                         all events. 
 * @param update_links : bool. Default false. If true update past_links sets
 *                      of all events.
 * @return Causet. The merged causet.  
 */
Causet Causet::merge(Causet pastCauset, Causet futCauset,
                           bool disjoint = false, bool update_futures = false,
                           bool update_links = false)
{
    int N = pastCauset.size();
    int M = futCauset.size();

    if (!disjoint)
    {
        for (set<int> p_f : futCauset._pasts)
        {
            // First update label of future events
            for (int i_f : p_f)
                {i_f += N;} 
            
            // Then add elements of pastCauset to past of future events
            // Do that with downwards loop, so that we always know that new
            // event goes at beginning of set ( O(n) instead of O(nlog(n)) )
            for (int i_p = N-1; i_p > -1; --i_p)
                {p_f.insert(p_f.begin(), i_p);}
        }
        // Then add elements of futCauset to future of past events if want
        // and past_links for future elements.
        if (update_futures && update_links)
        {
            for (set<int> pl_f : futCauset._past_links)
            {
                for (int i_p = N-1; i_p>-1; --i_p) //downwards for speed up
                {
                    set<int> f_p = pastCauset._futures[i_p];
                    // Add link between pastmost of fut and futmost of past
                    if (pl_f.empty() && f_p.empty())
                        {pl_f.insert(pl_f.begin(), i_p);}
                    // Add future event to future of past event
                    for (int i_f = 0; i_f < M; ++i_f)
                        {f_p.insert(f_p.end(), i_f);}
                }
            }
            vector<set<int>> pasts = {};
            pasts.insert(pasts.end(), pastCauset._pasts.begin(), 
                                      pastCauset._pasts.end());
            pasts.insert(pasts.end(), futCauset._pasts.begin(), 
                                      futCauset._pasts.end());
            vector<set<int>> futures = {};
            futures.insert(futures.end(), pastCauset._futures.begin(), 
                                      pastCauset._futures.end());
            futures.insert(futures.end(), futCauset._futures.begin(), 
                                      futCauset._futures.end());
            vector<set<int>> past_links = {};
            past_links.insert(past_links.end(), pastCauset._past_links.begin(), 
                                      pastCauset._past_links.end());
            past_links.insert(past_links.end(), futCauset._past_links.begin(), 
                                      futCauset._past_links.end());
            Causet Merged = Causet(pasts, futures, past_links);
            return Merged;
        }

        else if (update_futures)
        {
            for (set<int> f_p : pastCauset._futures)
            {
                // Add future event to future of past event
                for (int i_f = 0; i_f < M; ++i_f)
                    {f_p.insert(f_p.end(), i_f);}
            }
            vector<set<int>> pasts = {};
            pasts.insert(pasts.end(), pastCauset._pasts.begin(), 
                                      pastCauset._pasts.end());
            pasts.insert(pasts.end(), futCauset._pasts.begin(), 
                                      futCauset._pasts.end());
            vector<set<int>> futures = {};
            futures.insert(futures.end(), pastCauset._futures.begin(), 
                                      pastCauset._futures.end());
            futures.insert(futures.end(), futCauset._futures.begin(), 
                                      futCauset._futures.end());
            Causet Merged = Causet(pasts, futures, true);
            return Merged;
        }

        else if (update_links)
        {
            for (set<int> pl_f : futCauset._past_links)
            {
                for (int i_p = N-1; i_p>-1; --i_p) //downwards for speed up
                {
                    set<int> f_p = pastCauset._futures[i_p];
                    // Add link between pastmost of fut and futmost of past
                    if (pl_f.empty() && f_p.empty())
                        {pl_f.insert(pl_f.begin(), i_p);}
                }
            }
            vector<set<int>> pasts = {};
            pasts.insert(pasts.end(), pastCauset._pasts.begin(), 
                                      pastCauset._pasts.end());
            pasts.insert(pasts.end(), futCauset._pasts.begin(), 
                                      futCauset._pasts.end());
            vector<set<int>> past_links = {};
            past_links.insert(past_links.end(), pastCauset._past_links.begin(), 
                                      pastCauset._past_links.end());
            past_links.insert(past_links.end(), futCauset._past_links.begin(), 
                                      futCauset._past_links.end());
            Causet Merged = Causet(pasts, past_links, false);
            return Merged;
        }
        else
        {
            vector<set<int>> pasts = {};
            pasts.insert(pasts.end(), pastCauset._pasts.begin(), 
                                      pastCauset._pasts.end());
            pasts.insert(pasts.end(), futCauset._pasts.begin(), 
                                      futCauset._pasts.end());
            Causet Merged = Causet(pasts);
            return Merged;
        }
    } /*!disjoint*/
    else
    {
        // First update label of future events
        for (set<int> p_f : futCauset._pasts)
        {
            for (int i_f : p_f)
                {i_f += N;} 
        }
        // Then merge
        vector<set<int>> pasts = {};
        pasts.insert(pasts.end(), pastCauset._pasts.begin(), 
                                    pastCauset._pasts.end());
        pasts.insert(pasts.end(), futCauset._pasts.begin(), 
                                    futCauset._pasts.end());
        vector<set<int>> futures = {};
        futures.insert(futures.end(), pastCauset._futures.begin(), 
                                    pastCauset._futures.end());
        futures.insert(futures.end(), futCauset._futures.begin(), 
                                    futCauset._futures.end());
        vector<set<int>> past_links = {};
        past_links.insert(past_links.end(), pastCauset._past_links.begin(), 
                                    pastCauset._past_links.end());
        past_links.insert(past_links.end(), futCauset._past_links.begin(), 
                                    futCauset._past_links.end());
        Causet Merged = Causet(pasts, futures, past_links);
        return Merged;
    }
}







