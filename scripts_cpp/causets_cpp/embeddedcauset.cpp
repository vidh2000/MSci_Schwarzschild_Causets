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
 * @return : vector<bool> {x-y timelike, x<=y, x>y}.
 * @exception: returned if size of xvec and yvec diffrent than dimension of
 * spacetime.
 */
bool EmbeddedCauset::causality(vector<double> xvec, 
                                vector<double> yvec)
{
    auto xycausality = this->_spacetime.Causality();
    return xycausality(xvec, yvec, _spacetime._period, _spacetime._mass);
};


/**
 * @brief Causal relation according to the spacetime.
 * 
 * @param xvec: vector<double>, coordinates of x
 * @param yvec: vector<double>, coordinates of y
 * 
 * @return : vector<bool> {x-y timelike, x<=y, x>y}.
 * @exception: returned if size of xvec and yvec diffrent than dimension of
 * spacetime.
 */
std::vector<bool> EmbeddedCauset::general_causality(vector<double> xvec, 
                                                    vector<double> yvec)
{
    auto xycausality = this->_spacetime.General_Causality();
    return xycausality(xvec, yvec, _spacetime._period, _spacetime._mass);
};



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
bool EmbeddedCauset::areTimelike4D(vector<double> &xvec, vector<double> &yvec)
{
    double dt = (xvec[0]-yvec[0]);
    double dspacex = xvec[1]-yvec[1];
    double dspacey = xvec[2]-yvec[2];
    double dspacez = xvec[3]-yvec[3];
    // for(int i=1; i<dim; i++){
    //    dspace2 += (xvec[i]-yvec[i])*(xvec[i]-yvec[i]);
    // }
    if ((dt*dt)-(dspacex*dspacex + dspacey*dspacey + dspacez*dspacez)>0)
    {
        return true;}
    else{
        return false;}
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
    auto atimelikeb = _spacetime.General_Causality();
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

    _size -= labels.size();
    _dim = 0;
}


/**
 * @brief Get the interval between minimal and maximal element
 *          Essentially recreates the causet, getting rid of all
 *          the elements that aren't in the interval ->
 *          -> changes the pasts/futures sets and reduces the cmatrix
 * 
 *          The intervals are chosen at random and must contain #elemenets
 *          larger than min_size.
 * 
 * @param min_size - Minimal size of the interval (min # of elements in it) 
 * @param N_max - max number of tries to find the interval before stopping
 */
void EmbeddedCauset::get_interval(int min_size, int N_max) //=1000 max tries by default
{
    bool found = false; 
    int N_tries = 0;

    std::unordered_set<int> all_indices;
    for (int i = 0; i < _size; i++) {
        all_indices.insert(i);
  }

    while (!found)
    {
        // Failsafe
        if (N_tries > N_max){
            std::cout << "Couldn't find suitable interval in " << N_max
                << "tries" << std::endl;
            break;
        } 

        // Define mersenne_twister_engine Random Gen. (with random seed)
        std::random_device rd;
        int seed = rd();
        std::mt19937 gen(seed);
        std::uniform_real_distribution<> dis(0,_size);
        
        // Pick two random elements
        int e1 = (int) dis(gen), e2 =(int) dis(gen);
        int a; int b;
        if (e1 == e2){
            N_tries += 1;
            continue;
        }
        else if (e1 < e2){
            a = e1;
            b = e2;
        }
        else if (e1>e2){
            a = e2;
            b = e1;
        }
        else{
            N_tries += 1;
            continue;
        }
        int n = IntervalCard(a, b);
        if (n >= min_size && n<= _size)
        {  

            std::unordered_set<int> interval = set_intersection(
                        _futures[a], _pasts[b]);
            std::unordered_set<int> indices_to_remove =
                        set_diff(all_indices,interval);
            std::vector<int> to_discard(indices_to_remove.begin(),
                                        indices_to_remove.end());

            // Assume matrix is created and future and pasts but no links.
            EmbeddedCauset::discard(to_discard,true,true,false);
            found = true;
        }
        else{
            N_tries +=1;
            continue;
        }
    }
}



//////////////////////////////////////////////////////////////////////////////
//===========================================================================
// SAVE 
//===========================================================================
/////////////////////////////////////////////////////////////////////////////

/**
 * @brief Save causet attributes in file (ideally txt or csv)
 * =================================================================================
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




//////////////////////////////////////////////////////////////////////////////
//============================================================================
// MAKE ATTRIBUTES //=========================================================
//////////////////////////////////////////////////////////////////////////////
//============================================================================
/**
 * @brief Creates chosen attribues.Requires _size to have already be defined,
 *  and EVENTS ALREADY SORTED BY NATURAL LABELLING.
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
                std::cout<<"commented out\n";}
                //this->make_cmatrix_and_pastlinks(method, special);}
            else if (strcmp(sets_type, "future")==0){
                this->make_cmatrix_and_futlinks(method, special);}
        }
        else if (make_links == false && make_sets == true)
        {
            if (strcmp(sets_type, "past")==0){
                std::cout<<"commented out\n";}
                //this->make_cmatrix_and_pasts(method, special, use_transitivity);}
            else if (strcmp(sets_type, "future")==0){
                std::cout<<"commented out\n";}
                //this->make_cmatrix_and_futs(method, special, use_transitivity);}
        }
        else /*both make_sets and links*/
        {
            if (strcmp(sets_type, "past")==0){
                std::cout<<"commented out\n";}
                //this->make_cmatrix_and_allpasts(special);}
            else if (strcmp(sets_type, "future")==0){
                std::cout<<"commented out\n";}
                //this->make_cmatrix_and_allfuts(special);}
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
                std::cout<<"commented out\n";}
                //this->make_fut_links(method);}
        }
        else if (make_links == false && make_sets == true)
        {
            if (strcmp(sets_type, "past")==0){
                std::cout<<"commented out\n";}
                //this->make_pasts(method);}
            else if (strcmp(sets_type, "future")==0){
                std::cout<<"commented out\n";}
                //this->make_futures(method);}
        }
        else if (make_links == true && make_sets == true)
        {
            if (strcmp(sets_type, "past")==0){
                std::cout<<"commented out\n";}
                //this->make_all_pasts(method);}
            else if (strcmp(sets_type, "future")==0){
                std::cout<<"commented out\n";}
                //this->make_all_futures(method);}
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
            //std::cout << "Making cmatrix only with transitivity, without parallel inside\n";
            // std::cout << "j before i\n";
            // int special_factor = (special)? -1 : 1;
            // for(int j=1; j<_size; j++) //can skip the very first, i.e 0th
            // {
            //     for(int i=j-1; i>-1; i--) //i can only preceed j
            //     {
            //         if (_CMatrix[i][j] != 0){
            //             continue;}
            //         else
            //         {
            //             //if(xycausality(_coords[i],_coords[j],st_period,mass))
            //             {
            //                 _CMatrix[i][j] = special_factor;
            //                 // Obtain transitive relations
            //                 //#pragma omp parallel for //schedule(dynamic,8)
            //                 for (int k = i-1; k>-1; k--)
            //                 {
            //                     if(_CMatrix[k][i] != 0) //k<i<j -> k<j
            //                         { _CMatrix[k][j] = 1;}
            //                 }
            //             }
            //         }
            //     }
            // }
            //std::cout << "i before j\n";
            //int special_factor = (special)? -1 : 1;
            for(int i=1; i<_size; i++) //can skip the very first, i.e 0th
            {
                for(int j=i-1; j>-1; j--) //i can only preceed j
                {
                    if (_CMatrix[i][j] != 0){
                        continue;}
                    else
                    {
                        // note in flatspacetime funct, I set dim=4 hardcoded in spacetimes.cpp
                        if(xycausality(_coords[i],_coords[j],st_period,mass))
                        //if (areTimelike4D(_coords[i],_coords[j],4))
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
        else 
        {
            //std::cout << "Making cmatrix only without transitivity, without parallel\n";
            //std::cout << "j then i..\n";
            //#pragma omp parallel for// schedule(dynamic,8)
            //#pragma omp parallel for collapse(2)
            
            //std::cout << "with parallelisation\n";
            #pragma omp parallel for schedule(dynamic)
            for(int i=0; i<_size-1; i++) //can skip the very last, i.e Nth
            {
                for(int j=i+1; j<_size; j++) //i can only preceed j
                {
                    if(xycausality(_coords[i],_coords[j],st_period,mass))
                    //if(_spacetime.BH_causal4D(_coords[i],_coords[j],st_period,mass))
                    //if(_spacetime.Flat_causal(_coords[i],_coords[j],st_period,mass)[0])
                    //Returning a vector of boleans vs just boolean gives 33% boost in time
                    //if (areTimelike4D(_coords[i],_coords[j]))
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
 * - "Cmatrix": create from already existing _CMatrix
 * - "pasts": create from already existing futures
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
// void EmbeddedCauset::make_futures(const char* method)// = "coordinates")
// {   
//     _futures.resize(_size);
//     // Loop through coordinates t_min -> t_max
//     for (int i =_size-1; i>-1; i--)
//     {
//         //std::cout << "Event #"<< i+1 << std::endl;
//         for(int j=i+1; j<_size; j++)
//         {
//             if (set_contains(j, _futures[i]))
//                 {continue;}
//             else
//             {
//                 if (strcmp(method, "coordinates")==0 &&
//                     areTimelike(_coords[i], _coords[j]))
//                 {
//                     _futures[i].insert(j);
//                     _futures[i].insert(_futures[j].begin(), _futures[j].end());
//                 }
//             }
//         }
//     }
// }


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
            //std::cout << "Event #"<< i+1 << std::endl;
            for(int j=i+1; j<_size; j++) //j can only follow i
            {
                if (set_contains(j, _futures[i]))
                    {continue;}
                else
                {
                    if (xycausality(_coords[i],_coords[j],st_period,mass))
                        {_future_links[i].insert(j);}
                }
            }
        }
    }
    //std::cout <<"Finished sprinkling..." << std::endl;
}





///////////////////////////////////////////////////////////////////////////////
// Methods for counting links --> entropy
///////////////////////////////////////////////////////////////////////////////


/**
 * @brief First, creates from causal matrix _CMatrix kind of a set of 
 * future_links vectors such that if an element has one or more future links, 
 * then TWO ONLY, THE FIRST TWO, will be added to the vectors 
 * (as that is enough to see if an element is maximal).
 * Then counts links between maximal elements - below t_f and inside r_S - and 
 * maximal_but_one elements outside r_S. 

 * @param t_f Highest boundary for time. CURRENTLY UNUSED AS FIXED TO MAX.
 * @param r_S Schwarzschild radius

 * @return int Number of links
 */
int EmbeddedCauset::count_links_fromCMatrix(double& t_f, double r_S)
{
    if (strcmp(_spacetime._name, "BlackHole")==0)
    {
        if (_CMatrix.size()==0)
        {
            std::cout << "To create future link matrix, CMatrix must exist";
            throw std::invalid_argument("No CMatrix");}

        _future_links.resize(_size);
        
        #pragma omp parallel for schedule(dynamic)
        for (int i=0; i<_size; i++)
        {
            int n_links_of_i = 0; 
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
                            break;} //breaks k loop
                    }
                    if (!has_broken){
                        _future_links[i].insert(j);
                        n_links_of_i += 1;
                        if (n_links_of_i > 1)
                            {break;} /*breaks j loop, hence goes to next i*/
                    }
                }
            }
        }
        int N = this->count_links_BH(t_f,r_S);
        return N;
    }
    else /*Spacetime name not BlackHole*/
    {
        std::cout<<"Please choose 'BlackHole' for spacetime." <<
        "Other spacetimes might be available in the future."
        << std::endl;
        throw std::invalid_argument("Wrong spacetime");
    }
}




/**
 * @brief   Finds number of links in the causet connecting maximal elements 
 *          below t_f and inside the horizon with maximal-but-one elements
 *          outside the horizon.
 *          Currently works only for spacetime "Schwarzschild" in EForig coords, but
 *          could be expanded if needed in the future. 
 * 
 * @param t_f Highest boundary for time. CURRENTLY UNUSED AS FIXED TO MAX.
 * @param r_S Schwarzschild radius

 * @return int Number of links
 */
int EmbeddedCauset::count_links_BH(double& t_f, double r_S)
{
    if (!strcmp(_spacetime._name, "BlackHole")==0)
    {
        std::cout<<"Please choose 'BlackHole' for spacetime." <<
        "Other spacetimes might be available in the future" << std::endl;
        throw std::invalid_argument("Wrong spacetime");
    }

    // Number of links
    int N = 0;
    
    // To find point with lowest time component. Hypersurface set at t=0 btw.
    std::vector<double> min_times;
    std::vector<double> min_radii;
    std::vector<double> max_radii;
    for (int i = _size-2; i>-1; i--)
    {     
        if (_coords[i][1]>r_S) // i outside horizon
        {
            for (int j=i+1; j<_size; j++)
            {
                if (_coords[j][0]>_coords[i][0]) //t_j>t_i (SHOULD ALWAYS GO HERE)
                {
                    if (_coords[j][1]<r_S)  // j inside horizon
                    {
                        if (_future_links[j].size()==0 && // if j==maximal
                            _future_links[i].size()==1)   // if i links only to j  
                        {
                            // check if j-i is the link
                            if (set_contains(j,_future_links[i]))
                            {
                                N++;
                                min_times.push_back(_coords[i][0]);
                                min_radii.push_back(_coords[i][1]); //outside points
                                max_radii.push_back(_coords[j][1]); //inside points
                            }
                        }
                    }
                }
                else // t_j<t_i //
                {
                    std::cout << "t_j<t_i!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
                    if (_coords[i][1]<r_S && _coords[j][1]>r_S) // t_j<t_i
                    {
                        if (_future_links[i].size()==0 &&  // if i==maximal
                            _future_links[j].size()==1)  // if j links only to i  
                        {
                            // check if j-i is the link
                            if (set_contains(i,_future_links[j])) //faster if here
                            {
                                N++;
                                min_times.push_back(_coords[j][0]);
                            }
                        }
                    }
                }
            }
        }
    }
    double mintime = vecmin(min_times);
    double minradiusOut = vecmin(min_radii);
    double maxradiusOut = vecmax(min_radii);
    double maxradiusIn = vecmax(max_radii);
    double minradiusIn = vecmin(max_radii);
    std::cout << "t_min for elements in these links = " << mintime << std::endl; 
    std::cout << "r_max OUT = " << maxradiusOut << std::endl;
    std::cout << "r_min OUT = " << minradiusOut << std::endl;
    std::cout << "r_max IN = " << maxradiusIn << std::endl;
    std::cout << "r_min IN = " << minradiusIn << std::endl;
    return N;
}


/**
 * @brief First, creates from causal matrix _CMatrix KIND OF a set of 
 * future_links vectors such that if an element has one or more future links, 
 * then TWO ONLY, THE FIRST TWO, will be added to the vectors 
 * (as that is enough to see if an element is maximal).
 * Then counts lambdas between maximal 
 * elements -below t_f and inside r_S- and maximal_but_one elements outside r_S.
 * 
 * 
 * @param t_f Highest boundary for time. CURRENTLY UNUSED AS FIXED TO MAX.
 * @param r_S Schwarzschild radius

 * @return map<int, std::vector<int>> : key is maximal element label, value 
 *         is a vector of maximal but one connected elements' labels.
 */
std::map<int,std::vector<int>> EmbeddedCauset::get_lambdas(double& t_f, 
                                                            double r_S)
{
    if (strcmp(_spacetime._name, "BlackHole")==0)
    {
        
        if (_future_links.size() == _size) /*if already defined*/
        {
            return this->get_lambdas_from_futlinks(t_f,r_S);
        }

        else if (_CMatrix.size()==0)
        {
            std::cout << "To create future link matrix, CMatrix must exist";
            throw std::invalid_argument("No CMatrix");}
        else
        {
            _future_links.resize(_size);
        
            for (int i=0; i<_size; i++)
            {
                int n_links_of_i = 0;
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
                        if (!has_broken)
                        {
                            n_links_of_i += 1;
                            _future_links[i].insert(j);
                            if (n_links_of_i - 1 > 0)
                                {break;} /*breaks j loop, hence goes to next i*/
                        }
                    }
                }
            }
            return this->get_lambdas_from_futlinks(t_f,r_S);
        }
    }
    else /*Spacetime name not BlackHole*/
    {
        std::cout<<"Please choose 'BlackHole' for spacetime." <<
        "Other spacetimes might be available in the future."
        << std::endl;
        throw std::invalid_argument("Wrong spacetime");
    }
}


/**
 * @brief First, creates from causal matrix _CMatrix a kind of a set of 
 * future_links vectors such that if an element has one or more future links, 
 * then ONE AND ONLY ONE, THE FIRST, will be added to the vectors 
 * (as that is enough to see if an element is maximal).
 * Then counts lambdas between maximal elements -below t_f and inside r_S- 
 * and maximal_but_one elements outside r_S. 
 * 
 * @param t_f Highest boundary for time. CURRENTLY UNUSED AS FIXED TO MAX.
 * @param r_S Schwarzschild radius

 * @return map<int, int> : key is lambdas' size, value is number of
           such lambdas.
 *         Also, note, we also keep:
 *         result[-1] = mintime;
 *         result[-2] = innermost;
 *         result[-3] = outermost;
 */
std::map<int,double> EmbeddedCauset::count_lambdas(double& t_f, double r_S)
{
    if (strcmp(_spacetime._name, "BlackHole")==0)
    {
        if (_future_links.size() == _size) /*if already defined*/
        {
            std::cout<<"No need for futlinks, straight to counting molecules\n";
            std::map<int,double> sizes = get_lambdas_sizes(t_f,r_S);
            std::cout << "Finished get_lambdas_sizes" << std::endl;
            std::map<int,double> distr = get_lambdas_distr(sizes); 
            std::cout << "Finished get_lambdas_distr" << std::endl;
            return distr;
        }

        else if (_CMatrix.size()==0)
        {
            std::cout << "To create future link matrix, CMatrix must exist";
            throw std::invalid_argument("No CMatrix");
        }
        
        else
        {
            std::cout<<"Starting doing futlinks in count_lambdas"<<std::endl;
            _future_links.resize(_size);
        
            #pragma omp parallel
            for (int i=0; i<_size; i++)
            {
                int n_links_of_i = 0;
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
                        if (!has_broken)
                        {
                            n_links_of_i += 1;
                            _future_links[i].insert(j);
                            if (n_links_of_i - 1 > 0)
                                {break;} /*breaks j loop, hence goes to next i*/
                        }
                    }
                }
            }
            std::cout << "Finished done futlinks" << std::endl;
            std::map<int,double> sizes = get_lambdas_sizes(t_f,r_S);
            std::cout << "Finished get_lambdas_sizes" << std::endl;
            std::map<int,double> distr = get_lambdas_distr(sizes); 
            std::cout << "Finished get_lambdas_distr" << std::endl;
            return distr;
        }
    }
    else /*Spacetime name not BlackHole*/
    {
        std::cout<<"Please choose 'BlackHole' for spacetime." <<
        "Other spacetimes might be available in the future."
        << std::endl;
        throw std::invalid_argument("Wrong spacetime");
    }
}




/**
 * @brief First, creates from causal matrix _CMatrix KIND OF a set of 
 * future_links vectors such that if an element has 2 or more future links, 
 * then THREE ONLY, THE FIRST THREE, will be added to the vectors 
 * (as that is enough to see if an element is maximal but 2).
 * Then get Hawking's radiation molecules apb:
 * - p maximal but two below t_f and outside the horizon
 * - a maximal below t_f and outside the horizon
 * - b maximal below t_f and inside the hoizon
 * Note: 
 * - if a<b the molecule is "close";
 * - if b spacelike a the molecule is "open".
 * 
 * @param t_f Highest boundary for time. CURRENTLY UNUSED AS FIXED TO MAX.
 * @param r_S Schwarzschild radius

 * @return map<int, std::vector<int>> : key is minimal p element label, value 
           is the pair of labels of [a,b].
 */
std::map<int,std::vector<int>> EmbeddedCauset::get_HRVs(double& t_f, 
                                                        double r_S)
{
    if (strcmp(_spacetime._name, "BlackHole")==0)
    {
        if (_future_links.size() == _size) /*if already defined*/
        {
            return this->get_HRVs_from_futlinks(t_f,r_S);
        }
        else if (_CMatrix.size()==0)
        {
            std::cout << "To create future link matrix, CMatrix must exist";
            throw std::invalid_argument("No CMatrix");
        }
        else
        {
            _future_links.resize(_size);

            for (int i=0; i<_size; i++)
            {
                int n_links_of_i = 0;
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
                        if (!has_broken)
                        {
                            n_links_of_i += 1;
                            _future_links[i].insert(j);
                            if (n_links_of_i - 2 > 0)
                                {break;} /*breaks j loop, hence goes to next i*/
                        }
                    }
                }
            }
            return this->get_HRVs_from_futlinks(t_f,r_S);
        }
    }
    else /*Spacetime name not BlackHole*/
    {
        std::cout<<"Please choose 'BlackHole' for spacetime." <<
        "Other spacetimes might be available in the future."
        << std::endl;
        throw std::invalid_argument("Wrong spacetime");
    }
}


/**
 * @brief First, creates from causal matrix _CMatrix a kind of a set of 
 * future_links vectors such that if an element has one or more future links, 
 * then ONE AND ONLY ONE, THE FIRST, will be added to the vectors 
 * (as that is enough to see if an element is maximal).
 * Then counts lambdas between maximal elements -below t_f and inside r_S- 
 * and maximal_but_one elements outside r_S. 
 * 
 * @param t_f Highest boundary for time. CURRENTLY UNUSED AS FIXED TO MAX.
 * @param r_S Schwarzschild radius

 * @return map<int, int> : key is lambdas' size, value is number of
           such lambdas.
 *         Also, note, we also keep:
 *         result[-1] = mintime;
 *         result[-2] = innermost;
 *         result[-3] = outermost;
 */
std::map<int,double> EmbeddedCauset::count_HRVs(double& t_f, double r_S)
{
    if (strcmp(_spacetime._name, "BlackHole")==0)
    {
        if (_future_links.size() == _size) /*if already defined*/
        {
            return this->get_HRVs_distr_from_futlinks(t_f,r_S);
        }
        else if (_CMatrix.size()==0)
        {
            std::cout << "To create future link matrix, CMatrix must exist";
            throw std::invalid_argument("No CMatrix");
        }
        else
        {
            std::cout<<"Starting doing futlinks in count_HRVs"<<std::endl;
            _future_links.resize(_size);
        
            #pragma omp parallel
            for (int i=0; i<_size; i++)
            {
                int n_links_of_i = 0;
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
                        if (!has_broken)
                        {
                            n_links_of_i += 1;
                            _future_links[i].insert(j);
                            if (n_links_of_i - 2 > 0) //n links of i == 3
                                {break;} /*breaks j loop, hence goes to next i*/
                        }
                    }
                }
            }
            return this->get_HRVs_distr_from_futlinks(t_f,r_S);
        }
    }
    else /*Spacetime name not BlackHole*/
    {
        std::cout<<"Please choose 'BlackHole' for spacetime." <<
        "Other spacetimes might be available in the future."
        << std::endl;
        throw std::invalid_argument("Wrong spacetime");
    }
}


/**
 * @brief Save the following information in a file:
 * =================================================================================
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
 * ===================================== "sets" option    =======>
 * [6 to 6+size-1,:] pasts
 * [6+size+1 to 6+2size,:] futures
 * [6+2size+2 to 6+3size+1,:] past links
 * [6+3size+3 to 6+4size+2,:] future links
 * 
 * [6+size] -> "Coordinates"
 * [6+size+1 to 6+2size,:] -> coordinates
 * 
 * [6+2size+1] -> r_S
 * 
 * ==================================== "lambdas" option ========>
 * From [6+2size+2], for sets
 * "NLambdas_Sized_1", number of lambdas with size 1 (single links)
 * "NLambdas_Sized_2", number of lambdas with size 2 (the normal lambda)
 * "NLambdas_Sized_i", number of lambdas with size i (i-lambda)
 * etc for all sizes up to max...
 * 
 * From after that
 * "Lambda0" upvertex, downvertex1, downvertex2, etc... ENDL
 * "Lambda1", upvertex, downvertex2, downvertex2, etc... ENDL
 * etc for all lambdas...
 * 
 *  * ==================================== "HRVs" option ========>
 * From [6+2size+2], for sets
 * "NHRVs_Open", number of open HRVs (a spacelike b)
 * "NHRVs_Close", number of close HRVs (a<b)
 * 
 * From after that
 * "HRV0" p-vertex, upvertex_a, upvertex_b
 * "HRV1", p-vertex, upvertex_b, upvertex_b
 * etc for all HRVs...
 * 
 * 
 * @param path_file_ext 
 * @param storage_option const char*. Either "sets"(default) or "cmatrix".
 * @param t_f Highest boundary for time. CURRENTLY UNUSED AS FIXED TO MAX.
 * @param r_S Schwarzschild radius. Default 2.
 * @param molecule_option const char*. Either "lambdas"(default) or "HRVs".
 */
void EmbeddedCauset::save_molecules(const char* path_file_ext,
                                    const char* storage_option,
                                    double t_f, double r_S,
                                    const char* molecule_option)
{
    this->save_causet(path_file_ext, storage_option);

    std::fstream out;
    out.open(path_file_ext, std::ios::app);
    out<<std::endl<<"r_S," <<2*_spacetime._mass<<std::endl;
    _future_links.resize(0);

    if (strcmp(molecule_option, "lambdas")==0)
    {   
        int i = 0;
        std::cout<<"Getting the lambdas"<<std::endl;
        auto lambdas = get_lambdas(t_f, r_S);
        int N = lambdas.size();
        for (std::pair<int,std::vector<int>> lambda_i : lambdas)
        {
            out<<"Lambda"<<i<<",";
            out<<lambda_i.first<<",";
            for (int j = 0; j < lambda_i.second.size(); j++)
            {
                out<<lambda_i.second[j];
                if (j != lambda_i.second.size()-1)
                    out<<",";
            }
            if (i != N-1) out<<std::endl;
            i++;
        }
    }

    else if (strcmp(molecule_option, "HRVs")==0)
    {        
        int i = 0;
        auto HRVs = get_HRVs(t_f, r_S);
        int N = HRVs.size();
        for (std::pair<int,std::vector<int>> hrv_i : HRVs)
        {
            out<<"HRV"<<i<<",";
            out<<hrv_i.first<<",";
            for (int j = 0; j < hrv_i.second.size(); j++)
            {
                out<<hrv_i.second[j];
                if (j != hrv_i.second.size()-1)
                    out<<",";
            }
            if (i != N-1) out<<std::endl;
            i++;
        }
    }
    out.close();
}





// Counting Behind the Scenes

/**
 * @brief   Finds lambdas in the causet connecting maximal elements 
 *          below t_f and inside the horizon with maximal-but-one elements
 *          outside the horizon.
 *          Currently works only for spacetime "Schwarzschild" in EForig coords, but
 *          could be expanded if needed in the future. 
 * 
 * @param t_f Highest boundary for time. CURRENTLY UNUSED AS FIXED TO MAX.
 * @param r_S Schwarzschild radius

 * @return map<int, vector<int>> : key is label of maximal element, 
    value is vector of labels of maximal but one elements associated with it.
 */
std::map<int,std::vector<int>> EmbeddedCauset::get_lambdas_from_futlinks
                                                (double& t_f, double r_S)
{
    if (!strcmp(_spacetime._name, "BlackHole")==0)
    {
        std::cout<<"Please choose 'BlackHole' for spacetime." <<
        "Other spacetimes might be available in the future" << std::endl;
        throw std::invalid_argument("Wrong spacetime");
    }

    // Maps label of maximal element to size of its lambda
    std::map<int, std::vector<int>> lambdas;
    
    // To find point with lowest time component. Hypersurface set at t=tmax btw.
    double mintime = std::nan("");

    for (int j = 1; j<_size; j++)
    {
        // if j is maximal and inside the horizon
        if (_future_links[j].size()==0 && _coords[j][1]<r_S) 
        {
            lambdas[j] = {};
            for (int i = j-1; i>-1; i--)
            {
                if (_coords[j][0]>_coords[i][0]) //t_j>t_i SHOULD ALWAYS GO HERE
                {
                    // if i is maximal but one and outside the horizon
                    if (_future_links[i].size()==1 && _coords[i][1]>r_S)
                    {
                        if (set_contains(j,_future_links[i])) //i-j is link
                        {
                            lambdas[j].push_back(i);
                        }
                    }
                }
                else /* t_j<t_i */
                {
                    std::cout << "ERROR: t_j < t_i\n";
                }
            }
        }
    }
    return lambdas;
}


/**
 * @brief   Finds lambdas in the causet connecting maximal elements 
 *          below t_f and inside the horizon with maximal-but-one elements
 *          outside the horizon.
 *          Currently works only for spacetime "Schwarzschild" in EForig coords, but
 *          could be expanded if needed in the future. 
 * 
 * @param t_f Highest boundary for time. CURRENTLY UNUSED AS FIXED TO MAX.
 * @param r_S Schwarzschild radius

 * @return map<int, double> : key is label of maximal element, value is number of
           maximal but one elements associated with it.
 */
std::map<int,double> EmbeddedCauset::get_lambdas_sizes(double& t_f, double r_S)
{
    if (!strcmp(_spacetime._name, "BlackHole")==0)
    {
        std::cout<<"Please choose 'BlackHole' for spacetime." <<
        "Other spacetimes might be available in the future" << std::endl;
        throw std::invalid_argument("Wrong spacetime");
    }

    // Maps label of maximal element to size of its lambda
    std::map<int, double> lambdas;
    
    // To find point with lowest time component. Hypersurface set at t=tmax btw.
    std::vector<double> mintime_vec;
    std::vector<double> innermost_vec;
    std::vector<double> outermost_vec;
 
    for (int j = 1; j<_size; ++j)
    {
        // if j is maximal and inside the horizon
        if (_future_links[j].size()==0 && _coords[j][1]<r_S) 
        {
            lambdas[j] = 0;
            for (int i = j-1; i>-1; --i)
            {
                if (_coords[j][0]>_coords[i][0]) //t_j>t_i SHOULD ALWAYS GO HERE
                {
                    // if i is maximal but one and outside the horizon
                    if (_future_links[i].size()==1 && _coords[i][1]>r_S)
                    {
                        // if i-j is link
                        if (_future_links[i].find(j) != _future_links[i].end()) 
                        {
                            lambdas[j] += 1;
                            mintime_vec  .push_back(_coords[i][0]);
                            innermost_vec.push_back(_coords[j][1]);
                            outermost_vec.push_back(_coords[i][1]);
                        }
                    }
                }
                else /* t_j<t_i */
                {
                    std::cout << "ERROR: t_j < t_i\n";
                }
            }
        }
    }
    double mintime   = vecmin(mintime_vec);
    double innermost = vecmin(innermost_vec);
    double outermost = vecmax(outermost_vec);
    std::cout << "\nt_min for elements in lambdas = " << mintime << std::endl; 
    std::cout << "r_min for elements in lambdas = " << innermost<<std::endl; 
    std::cout << "r_max for elements in lambdas = " << outermost<<std::endl; 
    lambdas[-1] = mintime;
    lambdas[-2] = innermost;
    lambdas[-3] = outermost;
    return lambdas;
}


/**
 * @brief   Finds DISTRIBUTION of lambdas connecting maximal elements 
 *          below t_f and inside the horizon with maximal-but-one elements
 *          outside the horizon.
 *          Currently works only for spacetime "Schwarzschild" in EForig coords, 
 *          but could be expanded if needed in the future. 
 * 
 * @param lambdas map<int,double> : key is label of maximal element, value is 
 *        number of maximal but one elements associated with it, i.e. size.

 * @return map<int, double> : key is label of lambdas' size, value is number of
           occurrences.
 */
std::map<int,double> EmbeddedCauset::get_lambdas_distr(
    const std::map<int, double> &lambdas)
{
    // Maps label of maximal element to size of its lambda
    std::map<int, double> lambdas_distr;

    for (auto pair : lambdas)
    {
        if (pair.first > 0 && pair.second != 0)
        lambdas_distr[pair.second] += 1;

        else if (pair.first < 0)
        lambdas_distr[pair.first] = pair.second*1.;
    }
    return lambdas_distr;
}



/**
 * @brief   Finds HRVs in the causet.
 *          Currently works only for spacetime "Schwarzschild" in EForig coords, but
 *          could be expanded if needed in the future. 
 * 
 * @param t_f Highest boundary for time. CURRENTLY UNUSED AS FIXED TO MAX.
 * @param r_S Schwarzschild radius

 * @return map<int, vector<int>> : key is label of minimal element, 
    value is pair of labels of [a,b].
 */
std::map<int,std::vector<int>> EmbeddedCauset::get_HRVs_from_futlinks
                                                (double& t_f, double r_S)
{
    if (!strcmp(_spacetime._name, "BlackHole")==0)
    {
        std::cout<<"Please choose 'BlackHole' for spacetime." <<
        "Other spacetimes might be available in the future" << std::endl;
        throw std::invalid_argument("Wrong spacetime");
    }

    // Maps label of maximal element to size of its lambda
    std::map<int, std::vector<int>> HRVs;
    
    // To find point with lowest time component. Hypersurface set at t=tmax btw.
    std::vector<double> mintime_vec;
    std::vector<double> innermost_vec;
    std::vector<double> outermost_vec;
 
    for (int p = 1; p<_size; p++)
    {
        // if j is maximal and inside the horizon
        if (_future_links[p].size()==2 && _coords[p][1]>r_S) 
        {
            std::vector<int> ab;
            ab.insert(ab.end(),
                     _future_links[p].begin(), _future_links[p].end());
            int a = (ab[0] < ab[1])? ab[0] : ab[1];
            int b = (ab[0] < ab[1])? ab[1] : ab[0];
            if (_coords[b][1] < r_S && _coords[a][1] >= r_S)
            {
                //the one inside is maximal
                if (_future_links[b].size()==0) 
                {
                    //the one outside is either maximal 
                    // or maximal but one and only connected to b
                    if (_future_links[a].size()==0 ||
                        (_future_links[a].size()==1 && 
                         set_contains(b,_future_links[a]))
                        ) 
                    {
                        HRVs[p] = {a,b};
                        mintime_vec  .push_back(_coords[p][0]);
                        innermost_vec.push_back(_coords[b][1]);
                        outermost_vec.push_back(_coords[p][1]);
                        outermost_vec.push_back(_coords[a][1]);
                    }
                }
            }
        }
    }
    double mintime   = vecmin(mintime_vec);
    double innermost = vecmin(innermost_vec);
    double outermost = vecmax(outermost_vec);
    std::cout << "t_min for elements in the HRVs = " << mintime << std::endl; 
    std::cout << "r_min for elements in the HRVs = " << innermost<<std::endl; 
    std::cout << "r_max for elements in the HRVs = " << outermost<<std::endl; 
    return HRVs;
}


/**
 * @brief   Finds distribution of HRVs in the causet (how many open and close).
 *          Currently works only for spacetime "Schwarzschild" in EForig coords, but
 *          could be expanded if needed in the future. 
 * 
 * @param t_f Highest boundary for time. CURRENTLY UNUSED AS FIXED TO MAX.
 * @param r_S Schwarzschild radius

 * @return map<int, int> : key is 0 for open and 1 for close 
    value is pair of labels of [a,b]. Also, -3, -2, -1 for outermost, innermost
    and mintime.
 */
std::map<int,double> EmbeddedCauset::get_HRVs_distr_from_futlinks(double& t_f, 
                                                               double r_S)
{
    if (!strcmp(_spacetime._name, "BlackHole")==0)
    {
        std::cout<<"Please choose 'BlackHole' for spacetime." <<
        "Other spacetimes might be available in the future" << std::endl;
        throw std::invalid_argument("Wrong spacetime");
    }

    // Maps label of maximal element to size of its lambda
    std::map<int, double> HRVs_distr;
    
    // To find point with lowest time component. Hypersurface set at t=tmax btw.
    std::vector<double> mintime_vec;
    std::vector<double> innermost_vec;
    std::vector<double> outermost_vec;
 
    for (int p = 1; p<_size; p++)
    {
        // if j is maximal and inside the horizon
        if (_future_links[p].size()==2 && _coords[p][1]>r_S) 
        {
            std::vector<int> ab;
            ab.insert(ab.end(),
                     _future_links[p].begin(), _future_links[p].end());
            int a = (ab[0] < ab[1])? ab[0] : ab[1];
            int b = (ab[0] < ab[1])? ab[1] : ab[0];
            if (_coords[b][1] < r_S && _coords[a][1] >= r_S)
            {
                //the one inside is maximal
                if (_future_links[b].size()==0) 
                {
                    //the one outside is maximal -> open
                    if (_future_links[a].size()==0)
                    {
                        HRVs_distr[0] += 1;
                    }
                    // the one outside is only connected to a -> close
                    else if (_future_links[a].size()==1 && 
                             set_contains(b,_future_links[a]))
                    {
                        HRVs_distr[1] += 1;
                    }
                    
                    //update
                    mintime_vec  .push_back(_coords[p][0]);
                    innermost_vec.push_back(_coords[b][1]);
                    outermost_vec.push_back(_coords[p][1]);
                    outermost_vec.push_back(_coords[a][1]);
                }
            }
        }
    }
    double mintime   = vecmin(mintime_vec);
    double innermost = vecmin(innermost_vec);
    double outermost = vecmax(outermost_vec);
    std::cout << "t_min for elements in the HRVs = " << mintime << std::endl; 
    std::cout << "r_min for elements in the HRVs = " << innermost<<std::endl; 
    std::cout << "r_max for elements in the HRVs = " << outermost<<std::endl; 
    HRVs_distr[-1] = mintime;
    HRVs_distr[-2] = innermost;
    HRVs_distr[-3] = outermost;
    return HRVs_distr;
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