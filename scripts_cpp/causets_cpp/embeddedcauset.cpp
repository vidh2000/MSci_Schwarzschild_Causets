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
//EmbeddedCauset::EmbeddedCauset(){}

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
EmbeddedCauset::EmbeddedCauset(){}
EmbeddedCauset::EmbeddedCauset(Spacetime spacetime, 
                                CoordinateShape shape, 
                                vector<vector<double>> coordinates,
                                bool make_matrix,// = true,
                                bool special,// = false,
                                bool use_transitivity,// = true,
                                bool make_sets,// = false,
                                bool make_links,// = false,
                                const char* sets_type)// = "past
{
    _size = coordinates.size();
    _coords = coordinates;
    _spacetime = spacetime;
    _shape = shape;

    this->make_attrs("coordinates", make_matrix, special, use_transitivity,
                     make_sets, make_links, sets_type);
}




//=============================================================================
//=============================================================================
//MAKE ATTRS  //===============================================================
//=============================================================================
//=============================================================================

/**
 * @brief Creates chosen attribues.Requires _size to have already be defined,
 *  and events sorted by a possible natural labelling.
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
 * @param sets_type: const char* specifying the type of set:
 * - "past": (default) make past related sets 
 * - "future": make future related sets.
 */
void EmbeddedCauset::make_attrs (const char* method,// = "coordinates",
                                    bool make_matrix,
                                    bool special,// = false,
                                    bool use_transitivity,// = true,
                                    bool make_sets,// = false,
                                    bool make_links,// = false,
                                    const char* sets_type)// = "past")
{
    //std::cout << "Creating causet N=" << size << " via cmatrix" << std::endl;
    int special_factor = (special && (use_transitivity||make_links))? -1 : 1;
    _special_matrix = special && use_transitivity;

    if (make_matrix)
    {
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

    else /*Don't make matrix*/
    {
        if (make_links == true && make_sets == false)
        {
            if (strcmp(sets_type, "past")==0){
                this->make_past_links(method);}
            else if (strcmp(sets_type, "future")==0){
                this->make_fut_links(method);}
        }
        
        else if (make_links == false && make_sets == true)
        {
            if (strcmp(sets_type, "past")==0){
                this->make_pasts(method);}
            else if (strcmp(sets_type, "future")==0){
                this->make_futures(method);}
        }

        else if (make_links == true && make_sets == true)
        {
            if (strcmp(sets_type, "past")==0){
                this->make_all_pasts(method);}
            else if (strcmp(sets_type, "future")==0){
                this->make_all_futures(method);}
        }
        else
        {   
            std::cout<<"At least one among make_matrix, \
                        make_sets and make_links must be true"<<std::endl;
            throw std::invalid_argument("At least one among make_matrix, \
                                    make_sets and make_links must be true");
        }
    }
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


/**
 * @brief Compute Euclidean distances of points from _shape_center.
 * @return vector<double> : distances
 */
vector<double> EmbeddedCauset::eu_distances()
{
    vector<double> distances(_size, 0.0); 
    for (int i = 0; i < _size; i++)
    {
        vector<double>ivec = _coords[i];
        for (int mu = 0; mu < _spacetime._dim; mu++)
        {
            distances[i] += (ivec[mu] - _shape._center[mu])
                           *(ivec[mu] - _shape._center[mu]);
        }
        distances[i] = std::sqrt(distances[i]);
    }
    return distances;
}


/**
 * @brief Compute Euclidean distances of points from _shape_center.
 * @return vector<double> : distances
 */
vector<double> EmbeddedCauset::sp_radii()
{
    vector<double> radii(_size, 0.0); 
    for (int i = 0; i < _size; i++)
    {
        vector<double>ivec = _coords[i];
        for (int j = 1; j < _spacetime._dim; j++)
        {
            radii[i] += (ivec[j] - _shape._center[j])
                       *(ivec[j] - _shape._center[j]);
        }
        radii[i] = std::sqrt(radii[i]);
    }
    return radii;
}


/**
 * @brief Get maximum euclidean distance from center among sprinkled points.
 *
 * @return double : maximum value
 */
double EmbeddedCauset::max_eu_dist()
{
    double m = 0;
    for (int i = 0; i<_size; i++)
    {
        vector<double>ivec = _coords[i];
        double dist_i = 0;
        for (int mu = 0; mu < _spacetime._dim; mu++)
        {
            dist_i += (ivec[mu] - _shape._center[mu])
                     *(ivec[mu] - _shape._center[mu]);
        }
        if (dist_i>m) {m = dist_i*1.;}
    }
    return std::sqrt(m);
}


/**
 * @brief Get maximum SPATIAL distance from center among sprinkled points.
 *
 * @return double : maximum value
 */
double EmbeddedCauset::max_sp_rad()
{
    double m = 0;
    for (int i = 0; i<_size; i++)
    {
        vector<double>ivec = _coords[i];
        double rad_i = 0;
        for (int j = 1; j < _spacetime._dim; j++)
        {
            rad_i += (ivec[j] - _shape._center[j])
                     *(ivec[j] - _shape._center[j]);
        }
        if (rad_i>m) {m = rad_i*1.;}
    }
    return std::sqrt(m);
}


/**
 * @brief Get maximum value of coordinate "dim" among sprinkled points.
 * 
 * @param dim : int
 * @return double : maximum value
 */
double EmbeddedCauset::max_along(int dim)
{
    double m = _coords[0][dim];
    for (int i = 1; i<_size; i++)
    {
        if (_coords[i][dim]>m) {m = _coords[i][dim]*1.;}
    }
    return m;
}

/**
 * @brief Get minimum value of coordinate "dim" among sprinkled points.
 * 
 * @param dim : int
 * @return double : minimum value
 */
double EmbeddedCauset::min_along(int dim)
{
    double m = _coords[0][dim];
    for (int i = 1; i<_size; i++)
    {
        if (_coords[i][dim]<m) {m = _coords[i][dim]*1.;}
    }
    return m;
}




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
    //std::cout << "In areTimelike\n";
    auto atimelikeb = this->_spacetime.Causality();
    //std::cout << "After checking causality...\n";
    return atimelikeb(xvec, yvec, _spacetime._period, _spacetime._mass)[0];
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
    return atimelikeb(xvec, yvec, _spacetime._period, _spacetime._mass)[1];
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




//////////////////////////////////////////////////////////////////////////////
//============================================================================
// MAKE ATTRIBUTES BEHIND THE SCENES //=======================================
//////////////////////////////////////////////////////////////////////////////
//============================================================================

void EmbeddedCauset::make_cmatrix(const char* method,
                                    bool special,
                                    bool use_transitivity)
{
    if (strcmp(method, "coordinates")==0)
    {
        int special_factor = (special)? -1 : 1;
        _CMatrix.resize(_size, vector<int>(_size,0));
        std::cout << "Inside loop for making a cmatrix without links and sets..\n";
        if (use_transitivity)
        {
            for(int j=1; j<_size; j++) //can skip the very first, i.e 0th
            {
                for(int i=j-1; i>-1; i--) //i can only preceed j
                {
                    if (_CMatrix[i][j] != 0){
                        continue;}
                    else
                    {
                        if (this->areTimelike(_coords[i], _coords[j]))
                        {
                            _CMatrix[i][j] = special_factor;
                            for (int k = i-1; k>-1; k--)
                            {
                                if(_CMatrix[k][i] != 0) //k<i<j -> k<j
                                    { _CMatrix[k][j] = 1;}
                            }
                            
                        }
                    }
                }
            }
        }
        else 
        {
            for(int j=1; j<_size; j++) //can skip the very first, i.e 0th
            {
                for(int i=j-1; i>-1; i--) //i can only preceed j
                {
                    if (this->areTimelike(_coords[i], _coords[j])){
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
 * Transitivity is mandatory as links are bieng made. 
 */
void EmbeddedCauset::make_cmatrix_and_allpasts(bool special)
{
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
            if (areTimelike(_coords[i], _coords[j]))
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
 * Transitivity is mandatory as links are bieng made.  
 */
void EmbeddedCauset::make_cmatrix_and_allfuts(bool special)
{
    int special_factor = (special)? -1 : 1;
    _CMatrix.resize(_size, vector<int>(_size,0));
    _futures.resize(_size);
    _future_links.resize(_size);
    for(int i=_size-1; i>-1; i--) //can skip the very last
    {
        for(int j=i+1; j<_size; j++) //j can only follow i
        {
            if (_CMatrix[j][i] != 0)
                {continue;}
            if (areTimelike(_coords[i], _coords[j]))
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


void EmbeddedCauset::make_cmatrix_and_pasts(const char* method,
                                               bool special,
                                               bool use_transitivity)
{
    if (strcmp(method, "coordinates")==0)
    {
        int special_factor = (special && use_transitivity)? -1 : 1;
        _CMatrix.resize(_size, vector<int>(_size,0));
        _pasts.resize(_size);
        std::cout << "PASTS SIZE = " << _pasts.size() << std::endl; 
        for(int j=1; j<_size; j++) //can skip the very first, i.e 0th
        {
            for(int i=j-1; i>-1; i--) //i can only preceed j
            {
                if (_CMatrix[i][j] != 0){
                    continue;}
                else
                {
                    if (areTimelike(_coords[i], _coords[j]))
                    {
                        _CMatrix[i][j] = special_factor;
                        _pasts[j].insert(i);
                        _pasts[j].insert(_pasts[i].begin(), _pasts[i].end());
                        if (use_transitivity)
                        {
                            for (int k = i-1; k>-1; k--)
                            {
                                if(_CMatrix[k][i] != 0) //k<i<j -> k<j
                                    { _CMatrix[k][j] = 1;}
                            }
                        }
                    }
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


void EmbeddedCauset::make_cmatrix_and_futs(const char* method,
                                            bool special,
                                            bool use_transitivity)
{
    if (strcmp(method, "coordinates")==0)
    {
        int special_factor = (special && use_transitivity)? -1 : 1;
        _CMatrix.resize(_size, vector<int>(_size,0));
        _futures.resize(_size);
        for(int i=_size-1; i>-1; i--)
        {
            _CMatrix[i].resize(_size);
            for(int j=i+1; j<_size; j++)
            {
                if (_CMatrix[j][i] != 0)
                    {continue;}
                if (strcmp(method, "coordinates")==0 && 
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
    else
    {
        std::cout<<"Creation of Matrix and futures failed because currently\
        only method = 'coordinates' is supported."<<std::endl;
        throw std::invalid_argument("Only coordinates method currently \
        supported");
    }

}


void EmbeddedCauset::make_cmatrix_and_pastlinks(const char* method,
                                               bool special)
{
    if (strcmp(method, "coordinates")==0)
    {
        int special_factor = (special)? -1 : 1;
        _CMatrix.resize(_size, vector<int>(_size,0));
        _past_links.resize(_size);
        for(int j=1; j<_size; j++) //can skip the very first, i.e 0th
        {
            for(int i=j-1; i>-1; i--) //i can only preceed j
            {
                if (_CMatrix[i][j] != 0)
                    {continue;}
                if (areTimelike(_coords[i], _coords[j]))
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


void EmbeddedCauset::make_cmatrix_and_futlinks(const char* method,
                                               bool special)
{
    if (strcmp(method, "coordinates")==0)
    {
        int special_factor = (special)? -1 : 1;
        _CMatrix.resize(_size, vector<int>(_size,0));
        _future_links.resize(_size);
        for(int i=_size-1; i>-1; i--) //can skip the very last
        {
            for(int j=i+1; j<_size; j++) //j can only follow i
            {
                if (_CMatrix[j][i] != 0)
                    {continue;}
                if (areTimelike(_coords[i], _coords[j]))
                {
                    _CMatrix[i][j] = special_factor;
                    _future_links[i].insert(j);
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
    else
    {
        std::cout<<"Creation of Matrix and future_links failed as currently\
        only method = 'coordinates' is supported."<<std::endl;
        throw std::invalid_argument("Only coordinates method currently \
        supported");
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
                if (strcmp(method, "coordinates")==0 &&
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
                if (strcmp(method, "coordinates")==0 &&
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
                if (strcmp(method, "coordinates")==0 &&
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
                if (strcmp(method, "coordinates")==0 &&
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
                if (strcmp(method, "coordinates")==0 &&
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
                if (strcmp(method, "coordinates")==0 &&
                    areTimelike(_coords[i], _coords[j]))
                    {_future_links[i].insert(j);}
            }
        }
    }
    //std::cout <<"Finished sprinkling..." << std::endl;
}





// Destructor
EmbeddedCauset::~EmbeddedCauset(){}    

// Run with
// cd scripts_cpp/causet_cpp
// g++ -g causet.cpp shapes.cpp spacetimes.cpp embeddedcauset.cpp -std=c++17 -o embeddedcauset -O2
// .\embeddedcauset
// rm embeddedcauset.exe
// cd ../
// cd../
// int main(){
// std::cout << "embeddedcauset.cpp WORKS! :)";
// }