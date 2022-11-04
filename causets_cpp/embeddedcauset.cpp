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
//=============================================================================
//CONSTRUCTORS  //=============================================================
//=============================================================================
//=============================================================================
/**
 * @brief Creates embedded causet.
 * 
 * @param spacetime: Spacetime.
 * @param shape: CoordinateShape.
 * @param coordinates: vector<vector<double>>: ith entry is coordinates of 
 * event i.
 * @param make_matrix: bool, if true make matrix.
 * @param make_sets: bool, if true make set (see sets_type)
 * @param make_links: bool, if true make links (see sets_type)
 * @param sets_type: const char* specifying the type of set:
 * - "past": make _past_links
 * - "future": make _future_links
 * @param use_transitivity: bool, if true use also transitivity to establish
 * causality relations. 
 */
EmbeddedCauset::EmbeddedCauset(Spacetime spacetime, 
                                CoordinateShape shape, 
                                vector<vector<double>> coordinates,
                                bool make_matrix = true,
                                bool special = false,
                                bool use_transitivity = true,
                                bool make_sets = false,
                                bool make_links = false,
                                const char* sets_type = "past")
{
    _size = coordinates.size();
    _coords = coordinates;
    _spacetime = spacetime;
    _shape = shape;

    if (make_matrix)
        {
            this->make_cmatrix("coordinates", special, use_transitivity,
                                make_sets, make_links, sets_type);
        }

    else if (make_sets)
        {
            if (make_links)
            {
                if (sets_type == "past")
                    {this->make_all_pasts("coordinates");}
                else if (sets_type == "future")
                    {this->make_all_futures("coordinates");}
            }
            else
            {
                if (sets_type == "past")
                    {this->make_pasts("coordinates");}
                else if (sets_type == "future")
                    {this->make_futures("coordinates");}
            }
        }

    else if (make_links)
        {
            if (sets_type == "past")
                    {this->make_past_links("coordinates");}
            else if (sets_type == "future")
                {this->make_fut_links("coordinates");}   
        }
        
    else
        {
            throw std::invalid_argument("At least one among make_matrix, \
                                    make_sets and make_links must be true");
        }


}

//=============================================================================
//=============================================================================
//MAKE ATTRS  //===============================================================
//=============================================================================
//=============================================================================

/**
 * @brief Creates _pasts and _past_links, i.e. the sets of past and past 
 * links for each event. Requires _size to have already be defined, and events
 * sorted by a possible natural labelling.
 * 
 * @param method: const char*, possible choices are
 * - "coordinates": create from coordinates causality
 * - "pasts": create from already existing _pasts
 * - "futures": create from already existing _futures
 * @param special: bool, if true have C[i][j]=-1 if link
 * @param use_trasnitivity: bool, if true exploit transitivity
 * @param make_sets: bool, if true make set (see sets_type)
 * @param make_links: bool, if true make links (see sets_type)
 * @param sets_type: const char* specifying the type of set:
 * - "past": make _past_links
 * - "future": make _future_links
 */
void EmbeddedCauset::make_cmatrix (const char* method = "coordinates",
                                    bool special = false,
                                    bool use_transitivity = true,
                                    bool make_sets = false,
                                    bool make_links = false,
                                    const char* sets_type = "past")
{
    //std::cout << "Creating causet N=" << size << " via cmatrix" << std::endl;
    int8_t special_factor = (special)? -1 : 1;
    _special_matrix = special;
    _CMatrix.resize(_size);
    
    if (make_links == "none")
    {
        for(int i=0; i<_size; i++)
        {
            _CMatrix[i].resize(_size);
            for(int j=i-1; j>-1; j--)
            {
                if (_CMatrix[j][i] != 0)
                    {continue;}
                if (method == "coordinates" && areTimelike(_coords[i], _coords[j]))
                {
                    _CMatrix[j][i] = special_factor;
                    if (use_transitivity)
                    {
                        for (int k = j-1; k>-1; k--)
                        {
                            if(_CMatrix[k][j] != 0)
                                {_CMatrix[k][i] = 1;}
                        }
                    }
                }
            }
        }
    }

    else if (make_links == "past")
    {
        _past_links.resize(_size);
        for(int i=0; i<_size; i++)
        {
            _CMatrix[i].resize(_size);
            for(int j=i-1; j>-1; j--)
            {
                if (_CMatrix[j][i] != 0)
                    {continue;}
                if (method == "coordinates" && areTimelike(_coords[i], _coords[j]))
                {
                    _CMatrix[j][i] = special_factor;
                    _past_links[i].insert(j); 
                    if (use_transitivity)
                    {
                        for (int k = j-1; k>-1; k--)
                        {
                            if(_CMatrix[k][j] != 0)
                                {_CMatrix[k][i] = 1;}
                        }
                    }
                }
            }
        }
    }

    else if (make_links == "future")
    {
        _future_links.resize(_size);
        for(int i=_size-1; i>-1; i--)
        {
            _CMatrix[i].resize(_size);
            for(int j=i+1; j<_size; j++)
            {
                if (_CMatrix[j][i] != 0)
                    {continue;}
                if (method == "coordinates" && areTimelike(_coords[i], _coords[j]))
                {
                    _CMatrix[j][i] = special_factor;
                    _future_links[i].insert(j); 
                    if (use_transitivity)
                    {
                        for (int k = j+1; k<_size; k++)
                        {
                            if(_CMatrix[k][j] != 0)
                                {_CMatrix[k][i] = 1;}
                        }
                    }
                }
            }
        }
    }
}

/**
 * @brief Creates _pasts and _past_links, i.e. the sets of past and past 
 * links for each event. Requires _size to have already be defined, and events
 * sorted by a possible natural labelling.
 * 
 * @param method: const char*, possible choices are
 * - "coordinates": create from coordinates causality
 * - "Cmatrix": create from already existing _CMatrix
 * - "futures": create from already existing futures
 */
void EmbeddedCauset::make_all_pasts(const char* method = "coordinates")
{   
    _pasts.resize(_size);
    _past_links.resize(_size);
    // Loop through coordinates t_min -> t_max
    for (int i=1; i<_size; i++)
    {
        //std::cout << "Event #"<< i+1 << std::endl;
        for(int j=i-1; j>-1; j--)
        {
            // Check if j^th element is in pasts[i]
            if (set_contains(j, _pasts[i]))
                {continue;}
            else
            {
                if (method == "coordinates" &&
                    areTimelike(_coords[i], _coords[j]))
                {
                    _past_links[i].insert(j);
                    _pasts[i].insert(j);
                    _pasts[i].insert(_pasts[j].begin(), _pasts[j].end());
                }
            }
        }
    }
    //std::cout <<"Finished sprinkling..." << std::endl;
}

/**
 * @brief Creates _futures and _future_links, i.e. the sets of future and future 
 * links for each event. Requires _size to have already be defined, and events
 * sorted by a possible natural labelling.
 * 
 * @param method: const char*, possible choices are
 * - "coordinates": create from coordinates causality
 * - "Cmatrix": create from already existing _CMatrix
 * - "pasts": create from already existing futures
 */
void EmbeddedCauset::make_all_futures(const char* method = "coordinates")
{   
    _futures.resize(_size);
    _future_links.resize(_size);
    // Loop through coordinates t_min -> t_max
    for (int i =_size-1; i>-1; i--)
    {
        //std::cout << "Event #"<< i+1 << std::endl;
        for(int j=i+1; j<_size; j++)
        {
            if (set_contains(j, _futures[i]))
                {continue;}
            else
            {
                if (method == "coordinates" &&
                    areTimelike(_coords[i], _coords[j]))
                {
                    _future_links[i].insert(j);
                    _futures[i].insert(j);
                    _futures[i].insert(_futures[j].begin(), _futures[j].end());
                }
            }
        }
    }
    //std::cout <<"Finished sprinkling..." << std::endl;
}

/**
 * @brief Creates _pasts i.e. the sets of past events for each event. 
 * Requires _size to have already be defined, and events
 * sorted by a possible natural labelling.
 * 
 * @param method: const char*, possible choices are
 * - "coordinates": create from coordinates causality
 * - "Cmatrix": create from already existing _CMatrix
 * - "futures": create from already existing futures
 */
void EmbeddedCauset::make_pasts(const char* method = "coordinates")
{   
    _pasts.resize(_size);
    // Loop through coordinates t_min -> t_max
    for (int i=1; i<_size; i++)
    {
        //std::cout << "Event #"<< i+1 << std::endl;
        for(int j=i-1; j>-1; j--)
        {
            // Check if j^th element is in pasts[i]
            if (set_contains(j, _pasts[i]))
                {continue;}
            else
            {
                if (method == "coordinates" &&
                    areTimelike(_coords[i], _coords[j]))
                {
                    _pasts[i].insert(j);
                    _pasts[i].insert(_pasts[j].begin(), _pasts[j].end());
                }
            }
        }
    }
    //std::cout <<"Finished sprinkling..." << std::endl;
}

/**
 * @brief Creates _futures i.e. the sets of future events for each event. 
 * Requires _size to have already be defined, and events
 * sorted by a possible natural labelling.
 * 
 * @param method: const char*, possible choices are
 * - "coordinates": create from coordinates causality
 * - "Cmatrix": create from already existing _CMatrix
 * - "pasts": create from already existing pasts
 */
void EmbeddedCauset::make_futures(const char* method = "coordinates")
{   
    _futures.resize(_size);
    // Loop through coordinates t_min -> t_max
    for (int i =_size-1; i>-1; i--)
    {
        //std::cout << "Event #"<< i+1 << std::endl;
        for(int j=i+1; j<_size; j++)
        {
            if (set_contains(j, _futures[i]))
                {continue;}
            else
            {
                if (method == "coordinates" &&
                    areTimelike(_coords[i], _coords[j]))
                {
                    _futures[i].insert(j);
                    _futures[i].insert(_futures[j].begin(), _futures[j].end());
                }
            }
        }
    }
    //std::cout <<"Finished sprinkling..." << std::endl;
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
void EmbeddedCauset::make_past_links(const char* method = "coordinates")
{   
    _past_links.resize(_size);
    // Loop through coordinates t_min -> t_max
    for (int i=1; i<_size; i++)
    {
        //std::cout << "Event #"<< i+1 << std::endl;
        for(int j=i-1; j>-1; j--)
        {
            // Check if j^th element is in pasts[i]
            if (set_contains(j, _pasts[i]))
                {continue;}
            else
            {
                if (method == "coordinates" &&
                    areTimelike(_coords[i], _coords[j]))
                    {_past_links[i].insert(j);}
            }
        }
    }
    //std::cout <<"Finished sprinkling..." << std::endl;
}

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
void EmbeddedCauset::make_futures(const char* method = "coordinates")
{   
    _future_links.resize(_size);
    // Loop through coordinates t_min -> t_max
    for (int i =_size-1; i>-1; i--)
    {
        //std::cout << "Event #"<< i+1 << std::endl;
        for(int j=i+1; j<_size; j++)
        {
            if (set_contains(j, _futures[i]))
                {continue;}
            else
            {
                if (method == "coordinates" &&
                    areTimelike(_coords[i], _coords[j]))
                    {_future_links[i].insert(j);}
            }
        }
    }
    //std::cout <<"Finished sprinkling..." << std::endl;
}



//=============================================================================
//=============================================================================
//GETTERS     //===============================================================
//=============================================================================
//=============================================================================
int EmbeddedCauset::spacetime_dim()
    {return _spacetime._dim;}
double EmbeddedCauset:: density()
    {return _size/_shape.Volume();}
double EmbeddedCauset:: length_scale()
    {return std::pow( _size/_shape.Volume(), (double)1/_spacetime._dim );}



//=============================================================================
//=============================================================================
//RELATIONS   //===============================================================
//=============================================================================
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
//=============================================================================
//MODIFIERS   //===============================================================
//=============================================================================
//=============================================================================
void SprinkledCauset::add(std::set<CausetEvent> eventSet, bool unlink=false)
{
    double card_old = 1.0 * *this.Card();
    super().add(eventSet, unlink);  //??
    _intensity += (*this.Card() * 1.0 - card_old);
}


void SprinkledCauset::discard(std::set<CausetEvent> eventSet, bool unlink=false)
{
    double card_old = 1.0 * *this.Card();
    super().discard(eventSet, unlink);  //??
    _intensity *= (*this.Card() * 1.0 / card_old);
}