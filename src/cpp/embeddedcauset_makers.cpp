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
 * @brief Creates chosen attribues.Requires _size to have already be defined,
 *  and EVENTS ALREADY SORTED BY NATURAL LABELLING.
 * 
 * @note ONLY MAKE CMATRIX (make_matrix = true, everything else = false)
 *  is fastest. 
 * 
 * @param method: const char*, possible choices are
 * - "coordinates": create from coordinates causality
 * - "pasts": create from already existing _pasts
 * - "futures": create from already existing _futures
 * @param make_matrix : bool, if true(default) make _Cmatrix
 * @param special: bool, if true(default) have C[i][j]=-1 if link IFF also
 * use_transitivity is on.
 * @param use_transitivity: bool, if true(default) exploit transitivity. If 
 * make_links is true, it is compulsorily true.
 * @param make_sets: bool, if true (non-default) make set (see sets_type)
 * @param make_links: bool, if true (non-default) make links (see sets_type)
 * @param generation_mode: const char* specifying the generation method:
 * - "all with links" : make EVERYTHING overwriting make_matrix, make_sets,
 * make_links.
  * - "all" : Overwrite make_matrix and make_links and 
 * make CMatrix, _past and _futures. NOT LINKS. NOT SPECIAL. 
 * - "both only": makes both past and future, OVERWRITES make_matrix and
 * make_links, making them false. 
 * - "past": make _past_links
 * - "future": make _future_links
 */
void EmbeddedCauset::make_attrs (const char* method,// = "coordinates",
                                    bool make_matrix,
                                    bool special,// = false,
                                    bool use_transitivity,// = true,
                                    bool make_sets,// = false,
                                    bool make_links,// = false,
                                    const char* sets_type)// = "both only")
{
    // 1. Fix coordinates to EF(orginal) if it is BlackHole
    typedef void (*inversefunc)
    (std::vector<std::vector<double>>& coords, double mass, const char* EFtype);
    inversefunc inverse_transf = Spacetime::do_nothing;
    if (strcmp(_spacetime._name, "BlackHole")==0 
        && _spacetime._metricname!="EF(original)")
    {
        //change coordinates to EForig and save inverse function to go back
        inverse_transf = _spacetime.ToInEF_original(_coords);
        this->sort_coords(0, false);
    }

    // 2. Perform Causality
    if (strcmp(sets_type, "all with links")==0)
    {  
        this->make_cmatrix_and_allpasts(special);
        this->make_all_futures("coordinates");
    }
    
    else if (strcmp(sets_type, "all")==0)
    {
        this->make_all_but_links();
    }

    else if (strcmp(sets_type, "both only")==0)
    {
        this->make_sets(method);
    }

    else if (make_matrix)
    {
        _special_matrix = special && use_transitivity;
        if (make_links == false && make_sets == false)
        {
            this->make_cmatrix(method, special, use_transitivity);
        }
        else if (make_links == true && make_sets == false)
        {
            if (strcmp(sets_type, "past")==0){
                this->make_cmatrix_and_pastlinks(method, special);}
            else if (strcmp(sets_type, "future")==0){
                this->make_cmatrix_and_futlinks(method, special);}
        }
        else if (make_links == false && make_sets == true)
        {
            if (strcmp(sets_type, "past")==0){
                this->make_cmatrix_and_pasts(method, special, use_transitivity);}
            else if (strcmp(sets_type, "future")==0){
                this->make_cmatrix_and_futs(method, special, use_transitivity);}
        }
        else /*both make_sets and links*/
        {
            if (strcmp(sets_type, "past")==0){
                this->make_cmatrix_and_allpasts(special);}
            else if (strcmp(sets_type, "future")==0){
                this->make_cmatrix_and_allfuts(special);}
        }
    }

    else /*Don't make matrix neither both sets*/
    {
        if (make_links == true && make_sets == false)
        {
            if (strcmp(sets_type, "past")==0){
                std::cout<<"commented out\n";}
                //this->make_past_links(method);}
            else if (strcmp(sets_type, "future")==0){
                this->make_fut_links(method);}
        }
        else if (make_links == false && make_sets == true)
        {
            if (strcmp(sets_type, "past")==0){
                std::cout<<"commented out\n";}
                //this->make_pasts(method);}
            else if (strcmp(sets_type, "future")==0){
                this->make_futures(method);}
        }
        else if (make_links == true && make_sets == true)
        {
            if (strcmp(sets_type, "past")==0){
                std::cout<<"commented out\n";}
                //this->make_all_pasts(method);}
            else if (strcmp(sets_type, "future")==0){
                this->make_all_futures(method);}
        }
        else
        {   
            std::cout<<"Note: causet has no causal relations"<<std::endl;
        }
    }

    //Coords back to initial ones (does nothing if were not BlackHole or
    //were already EF(original))
    inverse_transf(_coords, _spacetime._mass, "original");
    return;
}


//////////////////////////////////////////////////////////////////////////////
//BEHIND THE SCENES

/**
 * @brief make CMatrix, pasts and futures from coordinates. NOT LINKS.
 * 
 */
void EmbeddedCauset::make_all_but_links()
{
    auto xycausality = this->_spacetime.Causality();
    std::vector<double> st_period = _spacetime._period;
    double mass = _spacetime._mass;

    _CMatrix.resize(_size, vector<int>(_size,0));
    _pasts.resize(_size);
    _futures.resize(_size);

    for(int j=1; j<_size; j++) 
    {
        for(int i=j-1; i>-1; i--) 
        {
            bool causalities = xycausality(_coords[i],_coords[j],
                                           st_period,mass);
            if (causalities) //i in past of j, j in future of i
            {
                _CMatrix[i][j] = 1;
                _pasts[j].insert(i);
                _pasts[j].insert(_pasts[i].begin(),_pasts[i].end());
                // Insert j into i's future and into
                // the future of elements in i's past
                _futures[i].insert(j);
                for (int ind_in_ipast : _pasts[i])
                {    
                    _futures[ind_in_ipast].insert(j);
                }
            }    
        }
    }
}


/**
 * @brief Makes _CMatrix from coordinates or past/futures set
 * 
 * @param method : const char* "coordinates" or "sets"
 * @param special : bool, if true, links identified with -1
 * @param use_transitivity : bool, if true, use transitivity to determine 
 * relations when possible.
 */
void EmbeddedCauset::make_cmatrix(const char* method,
                                    bool special,
                                    bool use_transitivity)
{ 
    auto xycausality = this->_spacetime.Causality();
    std::vector<double> st_period = _spacetime._period;
    double mass = _spacetime._mass;
    _CMatrix.resize(_size, vector<int>(_size,0));

    if (strcmp(method, "coordinates")==0)
    {
        if (use_transitivity)
        {
            for(int i=1; i<_size; i++) //can skip the very first, i.e 0th
            {
                for(int j=i-1; j>-1; j--) //i can only preceed j
                {
                    if (_CMatrix[i][j] != 0){
                        continue;}
                    else
                    {
                        if(xycausality(_coords[i],_coords[j],st_period,mass))
                        {
                            _CMatrix[i][j] = 1;//special_factor;
                            // Obtain transitive relations
                            //#pragma omp parallel for //schedule(dynamic,8)
                            for (int k = j-1; k>-1; k--)
                            {
                                if(_CMatrix[j][k] != 0) //k<i<j -> k<j
                                    { _CMatrix[i][k] = 1;}
                            }
                        }
                    }
                }
            }
        }
        else /*no transitivity*/
        {
            #pragma omp parallel for schedule(dynamic)
            for(int i=0; i<_size-1; i++) //can skip the very last, i.e Nth
            {
                for(int j=i+1; j<_size; j++) //i can only preceed j
                {
                    if(xycausality(_coords[i],_coords[j],st_period,mass))
                    {
                        _CMatrix[i][j] = 1;
                    }    
                }
            }
        }
    }
    else
    {
        std::cout<<"Creation of Matrix failed because currently\
        only method = 'coordinates' is supported."<<std::endl;
        throw std::invalid_argument("Only coordinates method currently \
        supported");
    }
}


/**
 * @brief Creates _CMatrix, _pasts and past_links. Can only be from coords. 
 * Transitivity is mandatory as links are being made. 
 */
void EmbeddedCauset::make_cmatrix_and_allpasts(bool special)
{
    auto xycausality = this->_spacetime.Causality();
    std::vector<double> st_period = _spacetime._period;
    double mass = _spacetime._mass;

    int special_factor = (special)? -1 : 1;
    _CMatrix.resize(_size, vector<int>(_size,0));
    _pasts.resize(_size);
    _past_links.resize(_size);
    for(int j=1; j<_size; j++) //can skip the very first, i.e 0th
    {
        for(int i=j-1; i>-1; i--) //i can only preceed j
        {
            if (_CMatrix[i][j] != 0)
                {continue;}
            if (xycausality(_coords[i],_coords[j],st_period,mass))
            {
                _CMatrix[i][j] = special_factor;
                _past_links[j].insert(i);
                _pasts[j].insert(i);
                _pasts[j].insert(_pasts[i].begin(), _pasts[i].end());
                // transitivity is mandatory if links are being made
                for (int k = i-1; k>-1; k--)
                {
                    if(_CMatrix[k][i] != 0) //k<i<j -> k<j
                        { _CMatrix[k][j] = 1;}
                }
            }
        }
    }
}


/**
 * @brief Creates _CMatrix, _futures and fut_links. Can only be from coords.
 * Transitivity is mandatory as links are being made.  
 */
void EmbeddedCauset::make_cmatrix_and_allfuts(bool special)
{
    auto xycausality = this->_spacetime.Causality();
    std::vector<double> st_period = _spacetime._period;
    double mass = _spacetime._mass;

    int special_factor = (special)? -1 : 1;
    if (_CMatrix.size()!=0) _CMatrix.clear();
    _CMatrix.resize(_size, vector<int>(_size,0));
    _futures.resize(_size);
    _future_links.resize(_size);
    for(int i=_size-2; i>-1; i--) //can skip the very last
    {
        for(int j=i+1; j<_size; j++) //j can only follow i
        {
            if (_CMatrix[i][j] != 0)
                {continue;}
            if (xycausality(_coords[i],_coords[j],st_period,mass))
            {
                _CMatrix[i][j] = special_factor;
                _future_links[i].insert(j);
                _futures[i].insert(j);
                _futures[i].insert(_futures[j].begin(), 
                                    _futures[j].end());
                // transitivity is mandatory if links are being made
                for (int k = j+1; k<_size; k++)
                {
                    if(_CMatrix[j][k] != 0) //i<j<k -> i<k
                        {_CMatrix[i][k] = 1;}
                }
            }
        }
    }
}


/**
 * @brief Make _CMatrix and _pasts
 * 
 * @param method const char* : either "coordinates" or "futures"
 * @param special bool : if true and use_transitivity, C_ij = -1 if ij link
 * @param use_transitivity bool : use transitivity where possible
 */
void EmbeddedCauset::make_cmatrix_and_pasts(const char* method,
                                               bool special,
                                               bool use_transitivity)
{
    if (strcmp(method, "coordinates")==0)
    {
        auto xycausality = this->_spacetime.Causality();
        std::vector<double> st_period = _spacetime._period;
        double mass = _spacetime._mass;

        _CMatrix.resize(_size, vector<int>(_size,0));
        _pasts.resize(_size); 
        if (use_transitivity)
        {
            int special_factor = (special)? -1 : 1;
            for(int j=1; j<_size; j++) //can skip the very first, i.e 0th
            {
                for(int i=j-1; i>-1; i--) //i can only preceed j
                {
                    if (_CMatrix[i][j] != 0){
                        continue;}
                    else if(xycausality(_coords[i],_coords[j],st_period,mass))
                    {
                        _CMatrix[i][j] = special_factor;
                        _pasts[j].insert(i);
                        _pasts[j].insert(_pasts[i].begin(), _pasts[i].end());
                        for (int k = i-1; k>-1; k--)
                        {
                            if(_CMatrix[k][i] != 0) //k<i<j -> k<j
                                { _CMatrix[k][j] = 1;}
                        }
                    }
                }
            }
        }
        else /*no transitivity*/
        for(int j=1; j<_size; j++) //can skip the very first, i.e 0th
        {
            for(int i=j-1; i>-1; i--) //i can only preceed j
            {
                if(xycausality(_coords[i],_coords[j],st_period,mass))
                {
                    _CMatrix[i][j] = 1;
                    _pasts[j].insert(i);
                }
            }
        }
    }
    else
    {
        std::cout<<"Creation of Matrix and pasts failed because currently\
        only method = 'coordinates' is supported."<<std::endl;
        throw std::invalid_argument("Only coordinates method currently \
        supported");
    }

}

/**
 * @brief Make _CMatrix and _futures
 * 
 * @param method const char* : either "coordinates" or "pasts"
 * @param special bool : if true and use_transitivity, C_ij = -1 if ij link
 * @param use_transitivity bool : use transitivity where possible
 */
void EmbeddedCauset::make_cmatrix_and_futs(const char* method,
                                            bool special,
                                            bool use_transitivity)
{
    if (strcmp(method, "coordinates")==0)
    {
        auto xycausality = this->_spacetime.Causality();
        std::vector<double> st_period = _spacetime._period;
        double mass = _spacetime._mass;

        _CMatrix.resize(_size, vector<int>(_size,0));
        _futures.resize(_size);
        if (use_transitivity)
        {
            int special_factor = (special && use_transitivity)? -1 : 1;
            for(int i=_size-1; i>-1; i--)
            {
                for(int j=i+1; j<_size; j++)
                {
                    if (_CMatrix[i][j] != 0)
                        {continue;}
                    else if(xycausality(_coords[i],_coords[j],st_period,mass))
                    {
                        _CMatrix[i][j] = special_factor;
                        _futures[i].insert(j);
                        _futures[i].insert(_futures[j].begin(), 
                                            _futures[j].end());
                        for (int k = j+1; k<_size; k++)
                        {
                            if(_CMatrix[k][j] != 0)
                                {_CMatrix[k][i] = 1;}
                        }
                    }
                }
            }
        }
        else /*no transitivity*/
        {
            for(int i=_size-1; i>-1; i--)
            {
                for(int j=i+1; j<_size; j++)
                {
                    if(xycausality(_coords[i],_coords[j],st_period,mass))
                    {
                        _CMatrix[i][j] = 1;
                        _futures[i].insert(j);
                    }
                }
            }
        }
    }
    else
    {
        std::cout<<"Creation of Matrix and futures failed because currently\
        only method = 'coordinates' is supported."<<std::endl;
        throw std::invalid_argument("Only coordinates method currently \
        supported");
    }
}


/**
 * @brief Make _CMatrix and _past_links. Transitivity is mandatory.
 * 
 * @param method const char* : either "coordinates" or "futures"
 * @param special bool : if true, C_ij = -1 if ij link
 */
void EmbeddedCauset::make_cmatrix_and_pastlinks(const char* method,
                                               bool special)
{
    if (strcmp(method, "coordinates")==0)
    {
        auto xycausality = this->_spacetime.Causality();
        std::vector<double> st_period = _spacetime._period;
        double mass = _spacetime._mass;

        int special_factor = (special)? -1 : 1;
        _CMatrix.resize(_size, vector<int>(_size,0));
        _past_links.resize(_size);

       
        for(int j=1; j<_size; j++) //can skip the very first, i.e 0th
        {
            
            for(int i=j-1; i>-1; i--) //i can only preceed j
            {
                if (_CMatrix[i][j] != 0)
                    {continue;}
                else if (xycausality(_coords[i],_coords[j],st_period,mass))
                {
                    _CMatrix[i][j] = special_factor;
                    _past_links[j].insert(i);
                    // transitivity is mandatory if links are being made
                    for (int k = i-1; k>-1; k--)
                    {
                        if(_CMatrix[k][i] != 0) //k<i<j -> k<j
                            { _CMatrix[k][j] = 1;}
                    }
                }
            }
        }
    }
    else
    {
        std::cout<<"Creation of Matrix and past_links failed because currently\
        only method = 'coordinates' is supported."<<std::endl;
        throw std::invalid_argument("Only coordinates method currently \
        supported");
    }

}


/**
 * @brief Make _CMatrix and _future_links.
 *        Transitivity is automatically applied, regardless of use_transitivity
 * 
 * @param method const char* : either "coordinates" or "pasts"
 * @param special bool : if true, C_ij = -1 if ij link
 */
void EmbeddedCauset::make_cmatrix_and_futlinks(const char* method,
                                               bool special)
{
    //std::cout << "Making cmatrix + futlinks, with parallel inside\n";
    if (strcmp(method, "coordinates")==0)
    {
        auto xycausality = this->_spacetime.Causality();
        std::vector<double> st_period = _spacetime._period;
        double mass = _spacetime._mass;

        int special_factor = (special)? -1 : 1;
        _CMatrix.resize(_size, vector<int>(_size,0));
        _future_links.resize(_size);

        for(int i=_size-1; i>-1; i--) //can skip the very last
        {
            //std::cout << "i="<<i<<"\n";
            for(int j=i+1; j<_size; j++) //j can only follow i
            {
                if (_CMatrix[i][j] != 0)
                    {continue;}
                else if (xycausality(_coords[i],_coords[j],st_period,mass))
                {
                    _CMatrix[i][j] = special_factor;
                    _future_links[i].insert(j);
                    // transitivity is mandatory if links are being made
                    //#pragma omp parallel for// schedule(dynamic,8)
                    for (int k = j+1; k<_size; k++)
                    {
                        if(_CMatrix[j][k] != 0) //i<j<k -> i<k
                            {_CMatrix[i][k] = 1;}
                    }
                }
            }           
        }
    }
    else
    {
        std::cout<<"Creation of Matrix and future_links failed as currently\
        only method = 'coordinates' is supported."<<std::endl;
        throw std::invalid_argument("Only coordinates method currently \
        supported");
    }
}


/**
 * @brief Make pasts and futures sets (not links).
 * 
 * @param method const char* : either "coordinates" or "cmatrix".
 */
void EmbeddedCauset::make_sets(const char* method)

{   
    _futures.resize(_size);
    _pasts.resize(_size);
    if (strcmp(method, "coordinates")==0)
    {
        // Define causality function pointer directly
        auto xycausality = this->_spacetime.Causality();
        std::vector<double> st_period = _spacetime._period;
        double mass = _spacetime._mass;

        // Loop through coordinates t_min -> t_max.
        // j>i automatically imposed as C_ij <-> i precedes j. 
        for (int j = 1; j<_size; j++)
        {
            for(int i=j-1; i>-1; i--)
            {

                // Does i precede j? == Is i<j?
                if (xycausality(_coords[i],_coords[j],st_period,mass)) 
                {
                    // Add i and its past to the past of j
                    _pasts[j].insert(i);
                    _pasts[j].insert(_pasts[i].begin(),_pasts[i].end());
                    // Insert j into i's future and into
                    // the future of elements in i's past
                    _futures[i].insert(j);
                    for (int ind_in_ipast : _pasts[i])
                    {    
                        _futures[ind_in_ipast].insert(j);
                    }
                }
            }
        }
    }
    //else if (strcmp(method, "cmatrix")==0)
    //    {this->make_sets_fromC();}
    else
    {
        std::cout<<"method must be 'coordinates' or 'cmatrix'"<<std::endl;
        throw std::invalid_argument("wrong method");
    }
}


/**
 * @brief Creates _pasts and _past_links, i.e. the sets of past and past 
 * links for each event. Requires _size to have already be defined, and events
 * sorted by a possible natural labelling.
 * 
 * @param method: const char*, possible choices are
 * - "coordinates": create from coordinates causality
 * - "cmatrix": create from already existing _CMatrix
 * - "futures": create from already existing futures
 */
// void EmbeddedCauset::make_all_pasts(const char* method)// = "coordinates")
// {   
//     _pasts.resize(_size);
//     _past_links.resize(_size);
//     // Loop through coordinates t_min -> t_max
//     for (int i=1; i<_size; i++)
//     {
//         //std::cout << "Event #"<< i+1 << std::endl;
//         for(int j=i-1; j>-1; j--)
//         {
//             // Check if j^th element is in pasts[i]
//             if (set_contains(j, _pasts[i]))
//                 {continue;}
//             else if (strcmp(method, "coordinates")==0 &&
//                     areTimelike(_coords[i], _coords[j]))
//             {
//                 _past_links[i].insert(j);
//                 _pasts[i].insert(j);
//                 _pasts[i].insert(_pasts[j].begin(), _pasts[j].end());
//             }
//         }
//     }
//     //std::cout <<"Finished sprinkling..." << std::endl;
// }


/**
 * @brief Creates _futures and _future_links, i.e. the sets of future and future 
 * links for each event. Requires _size to have already be defined, and events
 * sorted by a possible natural labelling.
 * 
 * @param method: const char*, possible choices are
 * - "coordinates": create from coordinates causality
 * OTHERS ARE NOT SUPPORTED
 * - "Cmatrix": create from already existing _CMatrix
 * - "pasts": create from already existing pasts
 */
void EmbeddedCauset::make_all_futures(const char* method)// = "coordinates")
{
    if (strcmp(method, "coordinates")==0)
    {
        auto xycausality = this->_spacetime.Causality();
        std::vector<double> st_period = _spacetime._period;
        double mass = _spacetime._mass;

        _futures.resize(_size);
        _future_links.resize(_size);
        for(int i=_size-2; i>-1; i--) //can skip the very last
        {
            for(int j=i+1; j<_size; j++) //j can only follow i
            {
                // Check if j^th element is in pasts[i]
                if (set_contains(j, _futures[i]))
                    {continue;}
                else if (xycausality(_coords[i],_coords[j],st_period,mass))
                {
                    _future_links[i].insert(j);
                    _futures[i].insert(j);
                    _futures[i].insert(_futures[j].begin(), 
                                        _futures[j].end());
                }
            }
        }
    }
    else
    {
        std::cout<<"method must be 'coordinates' or 'cmatrix', but, you know,"
                 <<" we have not implemented 'cmatrix'."<<std::endl;
        throw std::invalid_argument("wrong method");
    }
}



/**
 * @brief Creates _pasts i.e. the sets of past events for each event. 
 * Requires _size to have already be defined, and events
 * sorted by a possible natural labelling.
 * 
 * @param method: const char*, possible choices are
 * - "coordinates": create from coordinates causality
 * - "Cmatrix": create from already existing _CMatrix (not yet implemented)
 * - "futures": create from already existing futures (not yet implemented)
 */
// void EmbeddedCauset::make_pasts(const char* method)// = "coordinates")
// {   
//     _pasts.resize(_size);
//     // Loop through coordinates t_min -> t_max
//     for (int i=1; i<_size; i++)
//     {
//         //std::cout << "Event #"<< i+1 << std::endl;
//         for(int j=i-1; j>-1; j--)
//         {
//             // Check if j^th element is in pasts[i]
//             if (set_contains(j, _pasts[i]))
//                 {continue;}
//             else
//             {
//                 if (strcmp(method, "coordinates")==0 &&
//                     areTimelike(_coords[i], _coords[j]))
//                 {
//                     _pasts[i].insert(j);
//                     _pasts[i].insert(_pasts[j].begin(), _pasts[j].end());
//                 }
//             }
//         }
//     }
// }


/**
 * @brief Creates _futures i.e. the sets of future events for each event. 
 * Requires _size to have already be defined, and events
 * sorted by a possible natural labelling.
 * 
 * @param method: const char*, possible choices are
 * - "coordinates": create from coordinates causality
 * - "Cmatrix": create from already existing _CMatrix (not yet implemented)
 * - "pasts": create from already existing pasts (not yet implemented)
 */
void EmbeddedCauset::make_futures(const char* method)// = "coordinates")
{   
    auto xycausality = this->_spacetime.Causality();
    std::vector<double> st_period = _spacetime._period;
    double mass = _spacetime._mass;
    _futures.resize(_size);

    if (strcmp(method, "coordinates")==0)
    {
        #pragma omp parallel
        for(int i=0; i<_size-1; i++) //can skip the very last, i.e Nth
        {
            for(int j=i+1; j<_size; j++) //i can only preceed j
            {
                if(xycausality(_coords[i],_coords[j],st_period,mass))
                {
                    #pragma omp critical
                    _futures[i].insert(j);
                }    
            }
        }
    }
    else
    {
        std::cout<<"Creation of futures failed because currently\
        only method = 'coordinates' is supported."<<std::endl;
        throw std::invalid_argument("Only coordinates method currently \
        supported");
    }
}


/**
 * @brief Creates _past_links. 
 * Requires _size to have already be defined, and events
 * sorted by a possible natural labelling.
 * 
 * @param method: const char*, possible choices are
 * - "coordinates": create from coordinates causality
 * - "Cmatrix": create from already existing _CMatrix
 * - "futures": create from already existing futures
 */
// void EmbeddedCauset::make_past_links(const char* method)// = "coordinates")
// {   
//     _past_links.resize(_size);
//     // Loop through coordinates t_min -> t_max
//     for (int i=1; i<_size; i++)
//     {
//         //std::cout << "Event #"<< i+1 << std::endl;
//         for(int j=i-1; j>-1; j--)
//         {
//             // Check if j^th element is in pasts[i]
//             if (set_contains(j, _pasts[i]))
//                 {continue;}
//             else
//             {
//                 if (strcmp(method, "coordinates")==0 &&
//                     areTimelike(_coords[i], _coords[j]))
//                     {_past_links[i].insert(j);}
//             }
//         }
//     }
//     //std::cout <<"Finished sprinkling..." << std::endl;
// }


/**
 * @brief Creates _future_links. 
 * Requires _size to have already be defined, and events
 * sorted by a possible natural labelling.
 * 
 * @param method: const char*, possible choices are
 * - "coordinates": create from coordinates causality
 * - "Cmatrix": create from already existing _CMatrix
 * - "pasts": create from already existing pasts
 */
void EmbeddedCauset::make_fut_links(const char* method)// = "coordinates")
{   
    std::cout << "Making only futlinks, no parallel\n";
    _future_links.resize(_size);
    auto xycausality = this->_spacetime.Causality();
    std::vector<double> st_period = _spacetime._period;
    double mass = _spacetime._mass;
    if (strcmp(method, "coordinates")==0)
    {
        for (int i =_size-1; i>-1; i--) //can skip the very last
        {
            for(int j=i+1; j<_size; j++) //j can only follow i
            {
                if (xycausality(_coords[i],_coords[j],st_period,mass))
                    {_future_links[i].insert(j);}
            }
        }
    }
    else
    {
        std::cout<<"Creation of future links failed because currently\
        only method = 'coordinates' is supported."<<std::endl;
        throw std::invalid_argument("Only coordinates method currently \
        supported");
    }
}