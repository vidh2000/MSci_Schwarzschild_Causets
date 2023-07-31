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
 * @param special bool : if true, -1 if links.
 * @param use_transitivity: bool, if true use also transitivity to establish
 * causality relations. 
 * @param make_sets: bool, if true make set (see sets_type)
 * @param make_links: bool, if true make links (see sets_type)
 * @param generation_mode: const char* specifying the generation method:
 * - "all with links" : make EVERYTHING overwriting make_matrix, make_sets,
 * make_links.
 * - "both only": makes both past and future, OVERWRITES make_matrix and
 * make_links, making them false. 
 * - "past": make _past_links
 * - "future": make _future_links
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
                                const char* generation_mode)// = "both only"
{
    _size = coordinates.size();
    _coords = coordinates;
    _spacetime = spacetime;
    _shape = shape;

    this->make_attrs("coordinates", make_matrix, special, use_transitivity,
                     make_sets, make_links, generation_mode);
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


//////////////////////////////////////////////////////////////////////////////
//===========================================================================
// SAVE 
//===========================================================================
/////////////////////////////////////////////////////////////////////////////

/**
 * @brief Save causet attributes in file (ideally txt or csv)
 * ===========================================================================
 * [0,1] -> Storage Option
 * [1,1] -> size; 
 * [2,1] -> spacetime dimension;
 * [3,1] -> shape name;  
 * [4,1] -> spacetime name; 
 * 
 * ===================================== "cmatrix" option =======>
 * [5,0] -> "Matrix"
 * [6 to 6+size-1, 0 to size-1] -> Cmatrix; 
 * 
 * == "sets" option    =======>
 * [6 to 6+size-1,:] pasts
 * [6+size+1 to 6+2size,:] futures
 * [6+2size+2 to 6+3size+1,:] past links
 * [6+3size+3 to 6+4size+2,:] future links
 * 
 * [-size:] -> coordinates
 * 
 * @param path_file_ext const char* : path/file.ext 
 * @param storage_option const char* : cmatrix or sets
 */
void EmbeddedCauset::save_causet(const char* path_file_ext,
                                 const char* storage_option)
{
    std::ofstream out;
    out.open(path_file_ext);
    //if (!out.is_open()) std::cout<<"It is not open"<<std::endl;
    out<<"Storage option," << storage_option << std::endl;
    out<<"Size,"<<_size<<std::endl;
    out<<"Dimension,"<<_spacetime._dim<<std::endl;
    out<<"Shape,"<<_shape._name<<std::endl;
    out<<"Spacetime,"<<_spacetime._name<<std::endl;

    if (strcmp(storage_option, "cmatrix")==0)
    {
        out<<"Matrix,"<<std::endl;
        for (auto row : _CMatrix) 
        {
            for (auto col : row)
                {out << col <<',';}
            out<<std::endl;
        }
    }

    else if (strcmp(storage_option, "sets")==0)
    {
        if (!_pasts.size() || !_futures.size() ||
            !_past_links.size() || !_future_links.size())
            {
                std::cout << "You don't have all sets and links.\n";
                std::cout << "Please choose option 'all with links' \
                              for the causet generation" << std::endl;
                std::cout << "You're missing -> ";
                if (!_pasts.size()){
                    std::cout << "pasts" << std::endl;}
                if (!_futures.size()){
                    std::cout << "futures" << std::endl;}
                if (!_past_links.size()){
                    std::cout << "past links" << std::endl;}
                if (!_future_links.size()){
                    std::cout << "future links" << std::endl;}
                throw std::invalid_argument("dont have all sets+links");
            }
        out<<"Past sets, " << std::endl;
        for (auto past : _pasts)
        {
            for (auto e : past){
                out << e << ",";}
            out<<std::endl;
        }

        out<<"Future sets, " << std::endl;
        for (auto future : _futures)
        {
            for (auto e : future){
                out << e << ",";}
            out<<std::endl;
        }

        out<<"Past links sets, " << std::endl;
        for (auto past_links : _past_links)
        {
            for (auto e : past_links){
                out << e << ",";}
            out<<std::endl;
        }

        out<<"Future links sets, " << std::endl;
        for (auto fut_links : _future_links)
        {
            for (auto e : fut_links){
                out << e << ",";}
            out<<std::endl;
        }
    }

    else {
        std::cout << "Please choose 'cmatrix' or 'sets' option\n";
        throw std::invalid_argument("Choose right parameter");
    }

    out<<"Coordinates"<<std::endl;
    for (int i = 0; i < _size; i++) 
    {
        auto row = _coords[i];
        for (int mu = 0; mu < _spacetime._dim; mu++)
        {
            out<<row[mu];
            if (mu != _spacetime._dim -1)
                out<<",";
        }
        if (i != _size-1) 
            out<<std::endl;
    }

    out.close();
    return;
}