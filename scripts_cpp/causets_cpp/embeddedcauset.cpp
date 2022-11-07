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
#include "vecfunctions.h"
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
 * @brief Embed given coordinates in a causet.
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
                                bool make_matrix,// = true,
                                bool special,// = false,
                                bool use_transitivity,// = true,
                                bool make_sets,// = false,
                                bool make_links,// = false,
                                const char* sets_type)// = "past")
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
 * @param use_transitivity: bool, if true exploit transitivity. If make_links
 * is true, it is compulsorily true.
 * @param make_sets: bool, if true make set (see sets_type)
 * @param make_links: bool, if true make links (see sets_type)
 * @param sets_type: const char* specifying the type of set:
 * - "past": make _past_links
 * - "future": make _future_links
 */
void EmbeddedCauset::make_cmatrix (const char* method,// = "coordinates",
                                    bool special,// = false,
                                    bool use_transitivity,// = true,
                                    bool make_sets,// = false,
                                    bool make_links,// = false,
                                    const char* sets_type)// = "past")
{
    //std::cout << "Creating causet N=" << size << " via cmatrix" << std::endl;
    int8_t special_factor = (special)? -1 : 1;
    _special_matrix = special;
    _CMatrix.resize(_size, vector<int8_t>(_size));
    
    if (make_links == false && make_sets == false)
    {
        for(int i=0; i<_size; i++)
        {
            //_CMatrix[i].resize(_size);
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

    else if (make_links == true && make_sets == false)
    {
        if (sets_type == "past")
        {
            _past_links.resize(_size);
            for(int i=0; i<_size; i++)
            {
                for(int j=i-1; j>-1; j--)
                {
                    if (_CMatrix[j][i] != 0)
                        {continue;}
                    if (method == "coordinates" && areTimelike(_coords[i], _coords[j]))
                    {
                        _CMatrix[j][i] = special_factor;
                        _past_links[i].insert(j); 
                        //transitivity is mandatory if links are being done
                        for (int k = j-1; k>-1; k--)
                        {
                            if(_CMatrix[k][j] != 0)
                                {_CMatrix[k][i] = 1;}
                        }
                    }
                }
            }
        }
        else if (sets_type == "future")
        {
            _future_links.resize(_size);
            for(int i=_size-1; i>-1; i--)
            {
                for(int j=i+1; j<_size; j++)
                {
                    if (_CMatrix[j][i] != 0)
                        {continue;}
                    if (method == "coordinates" && areTimelike(_coords[i], _coords[j]))
                    {
                        _CMatrix[j][i] = special_factor;
                        _future_links[i].insert(j);
                        //transitivity is mandatory if links are being done
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
    
    else if (make_links == false && make_sets == true)
    {
        if (sets_type == "past")
        {
            _pasts.resize(_size);
            for(int i=0; i<_size; i++)
            {
                _CMatrix[i].resize(_size);
                for(int j=i-1; j>-1; j--)
                {
                    if (_CMatrix[j][i] != 0)
                        {continue;}
                    if (method == "coordinates" && 
                        areTimelike(_coords[i], _coords[j]))
                    {
                        _CMatrix[j][i] = special_factor;
                        _pasts[i].insert(j);
                        _pasts[i].insert(_pasts[j].begin(), _pasts[j].end());
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
        else if (sets_type == "future")
        {
            _futures.resize(_size);
            for(int i=_size-1; i>-1; i--)
            {
                _CMatrix[i].resize(_size);
                for(int j=i+1; j<_size; j++)
                {
                    if (_CMatrix[j][i] != 0)
                        {continue;}
                    if (method == "coordinates" && 
                        areTimelike(_coords[i], _coords[j]))
                    {
                        _CMatrix[j][i] = special_factor;
                        _futures[i].insert(j);
                        _futures[i].insert(_futures[j].begin(), 
                                           _futures[j].end());
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

    else /*both make_sets and links*/
    {
        if (sets_type == "past")
        {
            _pasts.resize(_size);
            for(int i=0; i<_size; i++)
            {
                _CMatrix[i].resize(_size);
                for(int j=i-1; j>-1; j--)
                {
                    if (_CMatrix[j][i] != 0)
                        {continue;}
                    if (method == "coordinates" && 
                        areTimelike(_coords[i], _coords[j]))
                    {
                        _CMatrix[j][i] = special_factor;
                        _past_links[i].insert(j);
                        _pasts[i].insert(j);
                        _pasts[i].insert(_pasts[j].begin(), _pasts[j].end());
                        // transitivity is mandatory if links are being made
                        for (int k = j-1; k>-1; k--)
                        {
                            if(_CMatrix[k][j] != 0)
                                {_CMatrix[k][i] = 1;}
                        }
                    }
                }
            }
        }
        else if (sets_type == "future")
        {
            _futures.resize(_size);
            for(int i=_size-1; i>-1; i--)
            {
                _CMatrix[i].resize(_size);
                for(int j=i+1; j<_size; j++)
                {
                    if (_CMatrix[j][i] != 0)
                        {continue;}
                    if (method == "coordinates" && 
                        areTimelike(_coords[i], _coords[j]))
                    {
                        _CMatrix[j][i] = special_factor;
                        _future_links[i].insert(j);
                        _futures[i].insert(j);
                        _futures[i].insert(_futures[j].begin(), 
                                           _futures[j].end());
                        // transitivity is mandatory if links are being made
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
void EmbeddedCauset::make_all_pasts(const char* method)// = "coordinates")
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
void EmbeddedCauset::make_all_futures(const char* method)// = "coordinates")
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
void EmbeddedCauset::make_pasts(const char* method)// = "coordinates")
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
void EmbeddedCauset::make_futures(const char* method)// = "coordinates")
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
void EmbeddedCauset::make_past_links(const char* method)// = "coordinates")
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
void EmbeddedCauset::make_fut_links(const char* method)// = "coordinates")
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
    _coords.erase(_coords.begin() + label);

    if (make_matrix)
    {
        if (_CMatrix.size())
        {
            _CMatrix.erase(_CMatrix.begin() + label);
            for (vector<int8_t> row : _CMatrix)
                {row.erase(row.begin() + label);}
        } 
    }
    if (make_sets)
    {
        if (_pasts.size())
        {
            _pasts.erase(_pasts.begin()+label);
            for (unordered_set<int> past_i : _pasts)
                {Causet::discard_from_set(past_i, label);}
        } 
        if (_futures.size())
        {
            _futures.erase(_futures.begin()+label);
            for (unordered_set<int> fut_i : _futures)
                {Causet::discard_from_set(fut_i, label);}
        }   
    }
    if (make_links)
    {
        if (_past_links.size())
        {
            _past_links.erase(_past_links.begin()+label);
            for (unordered_set<int> plinks_i : _past_links)
                {Causet::discard_from_set(plinks_i, label);}
        } 
        if (_future_links.size())
        {
            _future_links.erase(_future_links.begin()+label);
            for (unordered_set<int> flinks_i : _future_links)
                {Causet::discard_from_set(flinks_i, label);}
        }   
    }
    _size--;
    _dim = 0;
}


void EmbeddedCauset::discard(vector<int> labels,
                                    bool make_matrix, // = true, 
                                    bool make_sets, // = false,
                                    bool make_links) // = true)
{
    remove_indices(_coords, labels);

    if (make_matrix)
    {
        if (_CMatrix.size())
        {
            remove_indices(_CMatrix, labels);
            for (vector<int8_t> row : _CMatrix)
                {remove_indices(row, labels);}
        } 
    }
    if (make_sets)
    {
        if (_pasts.size())
        {
            remove_indices(_pasts, labels);
            for (unordered_set<int> past_i : _pasts)
                {Causet::discard_from_set(past_i, labels);}
        } 
        if (_futures.size())
        {
            remove_indices(_futures, labels);
            for (unordered_set<int> fut_i : _futures)
                {Causet::discard_from_set(fut_i, labels);}
        }   
    }
    if (make_links)
    {
        if (_past_links.size())
        {
            remove_indices(_past_links, labels);
            for (unordered_set<int> plinks_i : _past_links)
                {Causet::discard_from_set(plinks_i, labels);}
        } 
        if (_future_links.size())
        {
            remove_indices(_future_links, labels);
            for (unordered_set<int> flinks_i : _future_links)
                {Causet::discard_from_set(flinks_i, labels);}
        }   
    }
    _size--;
    _dim = 0;
}



int main(){
std::cout << "embeddedcauset.cpp WORKS! :)";

}