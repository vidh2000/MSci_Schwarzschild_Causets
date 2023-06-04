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
 * @brief If _futures are defined, then proceed with get_lambdas from futs.  
 * If _fut_links are defined, then proceed with get_lambdas_from_futlinks. 
 * If  _CMatrix is defined, creates from causal matrix KIND OF a set of 
 * futures vectors such that if an element has one or more futures, 
 * then TWO ONLY, THE FIRST TWO, will be added to the vectors 
 * (as that is enough to see if an element is maximal or maximal but 1).
 * 
 * Then gets lambdas between maximal elements -below t_f and inside r_S
 * - and maximal_but_one elements outside r_S.
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

        else if (_futures.size() == _size)
        {
            return this->get_lambdas_from_futs(t_f, r_S);
        }

        else if (_CMatrix.size()==0)
        {
            std::cout << "To create future link matrix, CMatrix must exist";
            throw std::invalid_argument("No CMatrix");}
        else
        {
            std::cout<<"Starting doing capped futs in count_lambdas"<<std::endl;
            _futures.resize(_size);
        
            #pragma omp parallel for
            for (int i=0; i<_size; i++)
            {
                int n_futs_of_i = 0;
                for (int j=i+1; j<_size; j++)
                {
                    #pragma omp critical
                    if (_CMatrix[i][j] == 1) {
                        _futures[i].insert(j);
                        n_futs_of_i += 1;
                    }

                    if (n_futs_of_i - 1 > 0){
                        break; /*break j loop, go to next i*/
                    }
                }
            }
            return this->get_lambdas_from_futs(t_f,r_S);
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
 * @brief If _futures are defined, then proceed with get_lambdas from futs.  
 * If _fut_links are defined, then proceed with get_lambdas_from_futlinks. 
 * If  _CMatrix is defined, creates from causal matrix KIND OF a set of 
 * futures vectors such that if an element has one or more futures, 
 * then TWO ONLY, THE FIRST TWO, will be added to the vectors 
 * (as that is enough to see if an element is maximal or maximal but 1).
 * 
 * Then gets lambdas between maximal elements -below t_f and inside r_S
 * - and maximal_but_one elements outside r_S.
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

        else if (_futures.size() == _size) /*if already defined*/
        {
            std::cout<<"No need for futs, straight to counting molecules\n";
            std::map<int,double> sizes = get_lambdas_sizes_from_futs(t_f,r_S);
            std::cout << "Finished get_lambdas_sizes" << std::endl;
            std::map<int,double> distr = get_lambdas_distr(sizes); 
            std::cout << "Finished get_lambdas_distr" << std::endl;
            return distr;
        }

        else if (_CMatrix.size()==0)
        {
            std::cout << "To create future, CMatrix must exist";
            throw std::invalid_argument("No CMatrix");
        }
        
        else
        {
            std::cout<<"Starting doing capped futs in count_lambdas"<<std::endl;
            _futures.resize(_size);
        
            //#pragma omp parallel for
            for (int i=0; i<_size; i++)
            {
                int n_futs_of_i = 0;
                for (int j=i+1; j<_size; j++)
                {
                    //#pragma omp critical
                    if (_CMatrix[i][j] == 1) {
                        _futures[i].insert(j);
                        n_futs_of_i += 1;
                    }

                    if (n_futs_of_i - 1 > 0){
                        break; /*break j loop, go to next i*/
                    }
                }
            }
            std::cout << "Finished done futs" << std::endl;
            std::map<int,double> sizes = get_lambdas_sizes_from_futs(t_f,r_S);
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
 * @brief If _futures are defined, then proceed with get_lambdas from futs.  
 * If  _CMatrix is defined, creates from causal matrix KIND OF a set of 
 * futures vectors such that if an element has one or more futures, 
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
        //  works wrong, so noooope
        // if (_future_links.size() == _size) /*if already defined*/
        // {
        //     return this->get_HRVs_from_futlinks(t_f,r_S);
        // }
        if (_futures.size() == _size){
            return this->get_HRVs_from_futs(t_f,r_S);
        }
        else if (_CMatrix.size()==0)
        {
            std::cout << "To create future link matrix, CMatrix must exist";
            throw std::invalid_argument("No CMatrix");
        }
        else
        {
            std::cout<<"Starting capped futures in count_HRVs"<<std::endl;
            _futures.resize(_size, {});
        
            //#pragma omp parallel for
            for (int i=0; i<_size-1; i++)
            {
                if (i==_size-2) std::cout<<"-"<<i<<"/"<<_size<<"-";
                int n_futs_of_i = 0;
                for (int j=i+1; j<_size; j++)
                {
                    //#pragma omp critical
                    if (_CMatrix[i][j] == 1) {
                        _futures[i].insert(j);
                        n_futs_of_i += 1;
                    }

                    if (n_futs_of_i - 2 > 0){ //n futs of i == 3
                        break; /*break j loop, go to next i*/
                    }
                }
            }
            return this->get_HRVs_from_futs(t_f,r_S);
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
 * @brief If _futures are defined, then proceed with get_lambdas from futs.  
 * If  _CMatrix is defined, creates from causal matrix KIND OF a set of 
 * futures vectors such that if an element has one or more futures, 
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

 * @return map<int, int> : key is HRV's type (0 open, 1 close), value is number of
           such molecules.
 *         Also, note, we also keep:
 *         result[-1] = mintime;
 *         result[-2] = innermost;
 *         result[-3] = outermost;
 */
std::map<int,double> EmbeddedCauset::count_HRVs(double& t_f, double r_S)
{
    if (strcmp(_spacetime._name, "BlackHole")==0)
    {
        if (_futures.size() == _size) /*if already defined*/
        {
            return this->get_HRVs_distr_from_futs(t_f,r_S);
        }
        else if (_CMatrix.size()==0)
        {
            std::cout << "To create future link matrix, CMatrix must exist";
            throw std::invalid_argument("No CMatrix");
        }
        else
        {
            // instead of futlinks, just do future, capped at 3 elements
            std::cout<<"Starting capped futures in count_HRVs"<<std::endl;
            _futures.resize(_size, {});
        
            //#pragma omp parallel for
            for (int i=0; i<_size-1; i++)
            {
                if (i==_size-2) std::cout<<"-"<<i<<"/"<<_size<<"-";
                int n_futs_of_i = 0;
                for (int j=i+1; j<_size; j++)
                {
                    //#pragma omp critical
                    if (_CMatrix[i][j] == 1) {
                        _futures[i].insert(j);
                        n_futs_of_i += 1;
                    }

                    if (n_futs_of_i - 2 > 0){ //n futs of i == 3
                        break; /*break j loop, go to next i*/
                    }
                }
            }
            std::cout<<"2314-";
            return this->get_HRVs_distr_from_futs(t_f,r_S);
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
 * [6+2size+1] -> "r_S", <r_S>
 * 
 * ==================================== "lambdas" option ========>
 * 
 * From after that
 * "Lambda0" upvertex, downvertex1, downvertex2, etc... ENDL
 * "Lambda1", upvertex, downvertex2, downvertex2, etc... ENDL
 * etc for all lambdas...
 * 
 * ==================================== "HRVs" option ========>
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




/**
 * @brief Save the following information in a file, given molecules
 * involving a total of N points were found. All indexes are scaled back to 
 * 0, 1, 2, ...
 * 
 * =================================================================================
 * [0,1] -> Storage Option
 * [1,1] -> size; 
 * [2,1] -> spacetime dimension;
 * [3,1] -> shape name;  
 * [4,1] -> spacetime name; 
 * 
 * [5,0] -> "Coordinates"
 * [6 to 6+N-1,:] -> coordinates
 * 
 * [6+N] -> "r_S", <r_S> 
 * 
 * ==================================== "lambdas" option ========>
 * From after that
 * "Lambda0" upvertex, downvertex1, downvertex2, etc... ENDL
 * "Lambda1", upvertex, downvertex2, downvertex2, etc... ENDL
 * etc for all lambdas...
 * 
 *  * ==================================== "HRVs" option ========>
 * From after that
 * "HRV0" p-vertex, upvertex_a, upvertex_b
 * "HRV1", p-vertex, upvertex_b, upvertex_b
 * etc for all HRVs...
 * 
 * 
 * @param path_file_ext 
 * @param t_f Highest boundary for time. CURRENTLY UNUSED AS FIXED TO MAX.
 * @param r_S Schwarzschild radius. Default 2.
 * @param molecule_option const char*. Either "lambdas"(default) or "HRVs".
 */
void EmbeddedCauset::save_molecules_only(const char* path_file_ext,
                                    double t_f, double r_S,
                                    const char* molecule_option)
{
    std::fstream out;
    out.open(path_file_ext);
    std::cout<<"Is open? in save_molecules_only"<<out.is_open()<<std::endl;
    //if (!out.is_open()) std::cout<<"It is not open"<<std::endl;
    out<<"Storage option," << "molecules only" << std::endl;
    out<<"Size,"<<_size<<std::endl;
    out<<"Dimension,"<<_spacetime._dim<<std::endl;
    out<<"Shape,"<<_shape._name<<std::endl;
    out<<"Spacetime,"<<_spacetime._name<<std::endl;

    if (strcmp(molecule_option, "lambdas")==0)
    {   
        int i = 0;
        std::cout<<"Getting the lambdas"<<std::endl;
        auto lambdas = get_lambdas(t_f, r_S);
        int N = lambdas.size();

        // Create and sort a vector of all points involved
        std::vector<int> all_points_involved;
        for (std::pair<int,std::vector<int>> lambda_i : lambdas)
        {
            all_points_involved.push_back(lambda_i.first);
            for (int j = 0; j < lambda_i.second.size(); j++)
            {
                all_points_involved.push_back(lambda_i.second[j]);
            }
            if (i != N-1) out<<std::endl;
            i++;
        }
        std::sort(all_points_involved.begin(), all_points_involved.end());

        // Create map from old label to new
        std::map<int, int> old_to_new_label;
        for (int i=0; i<all_points_involved.size(); i++){
            old_to_new_label[all_points_involved[i]] = i;
        }

        // Print coordinats
        out<<"Coordinates"<<std::endl;
        for (int n_label :  all_points_involved)
        {
            auto row = _coords[n_label];
            for (int mu = 0; mu < _spacetime._dim; mu++)
            {
                out<<row[mu];
                if (mu != _spacetime._dim -1)
                    out<<",";
            }
            if (i != _size-1) 
                out<<std::endl;
        }

        // Print lambdas in new labels
        out<<std::endl<<"r_S," <<2*_spacetime._mass<<std::endl;
        for (std::pair<int,std::vector<int>> lambda_i : lambdas)
        {
            out<<"Lambda"<<i<<",";
            out<<old_to_new_label[lambda_i.first]<<",";
            for (int j = 0; j < lambda_i.second.size(); j++)
            {
                out<<old_to_new_label[lambda_i.second[j]];
                if (j != lambda_i.second.size()-1)
                    out<<",";
            }
            if (i != N-1) out<<std::endl;
            i++;
        }
    }


    else if (strcmp(molecule_option, "HRVs")==0)
    {        
        std::cout<<"Doing HRVs in save_molecules_only"<<std::endl;
        int i = 0;
        auto HRVs = get_HRVs(t_f, r_S);
        int N = HRVs.size();

        // Create and sort a vector of all points involved
        std::vector<int> all_points_involved;
        for (std::pair<int,std::vector<int>> hrv_i : HRVs)
        {
            all_points_involved.push_back(hrv_i.first);
            for (int j = 0; j < hrv_i.second.size(); j++)
            {
                all_points_involved.push_back(hrv_i.second[j]);
            }
            if (i != N-1) out<<std::endl;
            i++;
        }
        std::sort(all_points_involved.begin(), all_points_involved.end());

        // Create map from old label to new
        std::map<int, int> old_to_new_label;
        for (int i=0; i<all_points_involved.size(); i++){
            old_to_new_label[all_points_involved[i]] = i;
        }

        // Print coordinats
        out<<"Coordinates"<<std::endl;
        for (int n_label :  all_points_involved)
        {
            auto row = _coords[n_label];
            for (int mu = 0; mu < _spacetime._dim; mu++)
            {
                out<<row[mu];
                if (mu != _spacetime._dim -1)
                    out<<",";
            }
            if (i != _size-1) 
                out<<std::endl;
        }

        // Print HRVs in new labels
        out<<std::endl<<"r_S," <<2*_spacetime._mass<<std::endl;
        for (std::pair<int,std::vector<int>> hrv_i : HRVs)
        {
            out<<"HRV"<<i<<",";
            out<<old_to_new_label[hrv_i.first]<<",";
            for (int j = 0; j < hrv_i.second.size(); j++)
            {
                out<<old_to_new_label[hrv_i.second[j]];
                if (j != hrv_i.second.size()-1)
                    out<<",";
            }
            if (i != N-1) out<<std::endl;
            i++;
        }
    }
    out.close();
}










////////////////////////////////////////////////////////////////////////
// Counting Behind the Scenes

/**
 * @brief   Finds lambdas in the causet connecting maximal elements 
 * below t_f and inside the horizon with maximal-but-one elements
 * outside the horizon.
 * Currently works only for spacetime "Schwarzschild" in EForig coords, but
 * could be expanded if needed in the future. 
 * 
 * @param t_f Highest boundary for time. CURRENTLY UNUSED AS FIXED TO MAX.
 * @param r_S Schwarzschild radius

 * @return map<int, vector<int>> : key is label of maximal element, 
    value is vector of labels of maximal but one elements associated with it.
 */
std::map<int,std::vector<int>> EmbeddedCauset::get_lambdas_from_futs
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

    for (int j = 1; j<_size; j++)
    {
        // if j is maximal and inside the horizon
        if (_futures[j].size()==0 && _coords[j][1]<r_S) 
        {
            lambdas[j] = {};
            for (int i = j-1; i>-1; i--)
            {
                if (_coords[j][0]>_coords[i][0]) //t_j>t_i SHOULD ALWAYS GO HERE
                {
                    // if i is maximal but one and outside the horizon
                    if (_futures[i].size()==1 && _coords[i][1]>r_S)
                    {
                         //i-j is link
                        if (_futures[i].find(j) != _futures[i].end())
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
std::map<int,double> EmbeddedCauset::get_lambdas_sizes_from_futs(double& t_f, 
                                                                double r_S)
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
 
    //#pragma omp parallel for //schedule(dynamic)
    for (int j = 1; j<_size; ++j)
    {
        // if j is maximal and inside the horizon
        if (_futures[j].size()==0 && _coords[j][1]<r_S) 
        {
            lambdas[j] = 0;
            for (int i = j-1; i>-1; --i)
            {
                // if i is maximal but one and outside the horizon
                if (_futures[i].size()==1 && _coords[i][1]>r_S)
                {
                    // if i-j is link
                    if (_futures[i].find(j) != _futures[i].end()) 
                    {
                        //#pragma omp critical
                        {
                        lambdas[j] += 1;
                        mintime_vec  .push_back(_coords[i][0]);
                        innermost_vec.push_back(_coords[j][1]);
                        outermost_vec.push_back(_coords[i][1]);
                        }
                    }
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
 * @brief  Finds HRVs in the causet, given futures of points.
 * Currently works only for spacetime "Schwarzschild" in EForig coords, but
 * could be expanded if needed in the future. 
 * 
 * @param t_f Highest boundary for time. CURRENTLY UNUSED AS FIXED TO MAX.
 * @param r_S Schwarzschild radius

 * @return map<int, vector<int>> : key is label of minimal element, 
    value is pair of labels of [a,b].
 */
std::map<int,std::vector<int>> EmbeddedCauset::get_HRVs_from_futs
                                                (double& t_f, double r_S)
{
    if (!strcmp(_spacetime._name, "BlackHole")==0)
    {
        std::cout<<"Please choose 'BlackHole' for spacetime." <<
        "Other spacetimes might be available in the future" << std::endl;
        throw std::invalid_argument("Wrong spacetime");
    }

    // Maps label of minimal element to top elements
    std::map<int, std::vector<int>> HRVs;
    
    for (int p = 0; p<_size; p++)
    {
        // if p is maximal but 2 outside the horizon
        if (_futures[p].size()==2 && _coords[p][1]>r_S) 
        {
            // p<a<b (number wise, not necessarily causality wise for a and b)
            std::vector<int> ab;
            ab.insert(ab.end(),
                     _futures[p].begin(), _futures[p].end());
            int a = (ab[0] < ab[1])? ab[0] : ab[1];
            int b = (ab[0] < ab[1])? ab[1] : ab[0];

            HRVs[p] = {a,b};
        }
    }
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
std::map<int,double> EmbeddedCauset::get_HRVs_distr_from_futs(double& t_f, 
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
 
    for (int p = 0; p<_size; p++)
    {
        // if p is maximal but 2 outside the horizon
        if (_futures[p].size()==2 && _coords[p][1]>r_S) 
        {
            // p<a<b (number wise, not necessarily causality wise for a and b)
            std::vector<int> ab;
            ab.insert(ab.end(),
                     _futures[p].begin(), _futures[p].end());
            int a = (ab[0] < ab[1])? ab[0] : ab[1];
            int b = (ab[0] < ab[1])? ab[1] : ab[0];

            // b in and a out
            if (_coords[b][1] <= r_S && _coords[a][1] > r_S)
            {
                // the one inside is authomatically maximal
                // if a is maximal
                // -> open HRV
                if (_futures[a].size()==0)
                {
                    HRVs_distr[0] += 1;
                    mintime_vec  .push_back(_coords[p][0]);
                    innermost_vec.push_back(_coords[b][1]);
                    outermost_vec.push_back(_coords[p][1]);
                    outermost_vec.push_back(_coords[a][1]);
                }
                // if a not maximal, it can only be connected to b
                // -> close HRV 
                else
                {
                    HRVs_distr[1] += 1;
                    mintime_vec  .push_back(_coords[p][0]);
                    innermost_vec.push_back(_coords[b][1]);
                    outermost_vec.push_back(_coords[p][1]);
                    outermost_vec.push_back(_coords[a][1]);
                }
            }

            // b out and a in
            if (_coords[a][1] <= r_S && _coords[b][1] > r_S)
            {
                //the one inside is maximal (authomatically implied)
                // b is larger and out, so can only be open
                HRVs_distr[0] += 1;
                mintime_vec  .push_back(_coords[p][0]);
                innermost_vec.push_back(_coords[a][1]);
                outermost_vec.push_back(_coords[p][1]);
                outermost_vec.push_back(_coords[b][1]);
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
 * @brief  Finds distribution of HRVs in the causet (how many open and close).
 * Currently works only for spacetime "Schwarzschild" in EForig coords, but
 * could be expanded if needed in the future. 
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

    // Maps type of HRV to its number
    std::map<int, double> HRVs_distr;
    
    // To find point with lowest time component. Hypersurface set at t=tmax btw.
    std::vector<double> mintime_vec;
    std::vector<double> innermost_vec;
    std::vector<double> outermost_vec;
 
    for (int p = 0; p<_size-2; p++)
    {
        // if p has exactly 2 links: these are a and b, one in and one out.
        // If both are maximal we have open HRV 
        if (_future_links[p].size()==2 && _coords[p][1]>r_S) 
        {
            // set a<b
            std::vector<int> ab;
            ab.insert(ab.end(),
                     _future_links[p].begin(), _future_links[p].end());
            int a = (ab[0] < ab[1])? ab[0] : ab[1];
            int b = (ab[0] < ab[1])? ab[1] : ab[0];

            // both a and b must be maximal
            if (_future_links[b].size()==0 && _future_links[a].size()==0) 
            {
                // if b in and a is out, we got our opn HRV
                if ((_coords[b][1]<=r_S && _coords[a][1]>r_S))
                {
                    HRVs_distr[0] += 1;
                    mintime_vec  .push_back(_coords[p][0]);
                    innermost_vec.push_back(_coords[b][1]);
                    outermost_vec.push_back(_coords[p][1]);
                    outermost_vec.push_back(_coords[a][1]);
                }
                // if b out and a in, we got our opn HRV
                else if (_coords[a][1]<=r_S && _coords[b][1]>r_S)
                {
                    HRVs_distr[0] += 1;
                    mintime_vec  .push_back(_coords[p][0]);
                    innermost_vec.push_back(_coords[a][1]);
                    outermost_vec.push_back(_coords[p][1]);
                    outermost_vec.push_back(_coords[b][1]);
                }
            }
        }

        // if p has exactly 1 link, this is a
        // If a is out and has exactly one link to b, maximal inside,
        // we have close HRV 
        else if (_future_links[p].size()==1 && _coords[p][1]>r_S) 
        {
            // get a
            std::vector<int> avec;
            avec.insert(avec.end(),
                     _future_links[p].begin(), _future_links[p].end());
            int a = avec[0];

            // if a is out and has exactly one link
            if ( _coords[a][1]>r_S && _future_links[a].size()==1)
            {
                // get b
                std::vector<int> bvec;
                avec.insert(bvec.end(),
                        _future_links[a].begin(), _future_links[a].end());
                int b = bvec[0];

                // if b is in and maximal we got close HRV
                if (_coords[b][1]<=r_S && _future_links[b].size()==0) 
                {
                    HRVs_distr[1] += 1;
                    mintime_vec  .push_back(_coords[p][0]);
                    innermost_vec.push_back(_coords[b][1]);
                    outermost_vec.push_back(_coords[p][1]);
                    outermost_vec.push_back(_coords[a][1]);
                }
            }
        }
    }
    double mintime   = 0;
    double innermost = 0;
    double outermost = 0;
    if (mintime_vec.size()>0){
    auto mintime_iter   = std::min_element(mintime_vec.begin(), 
                                           mintime_vec.end());
    mintime = *mintime_iter;
    }
    if (innermost_vec.size()>0){
    auto innermost_iter = std::min_element(innermost_vec.begin(), 
                                         innermost_vec.end());
    innermost = *innermost_iter;
    }
    if (outermost_vec.size()>0){
    auto outermost_iter = std::max_element(outermost_vec.begin(), 
                                         outermost_vec.end());
    outermost = *outermost_iter;
    }
    std::cout << "t_min for elements in the HRVs = " << mintime << std::endl; 
    std::cout << "r_min for elements in the HRVs = " << innermost<<std::endl; 
    std::cout << "r_max for elements in the HRVs = " << outermost<<std::endl; 
    HRVs_distr[-1] = mintime;
    HRVs_distr[-2] = innermost;
    HRVs_distr[-3] = outermost;
    return HRVs_distr;
}




/**
 * @brief  IS WRONG!!!!!!!!!!!!!!!!!! Finds HRVs in the causet.
 * Currently works only for spacetime "Schwarzschild" in EForig coords, but
 * could be expanded if needed in the future. 
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
                         _future_links[a].count(a)==1)
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
 * @brief First, creates from causal matrix _CMatrix a vector of capped
 * futures sets such that if an element has one or more future elements, 
 * then ONE AND ONLY ONE, THE FIRST, will be added to the vectors 
 * (as that is enough to see if an element is maximal).
 * Then counts lambdas between maximal elements -below t_f and inside r_S- 
 * and maximal_but_one elements outside r_S. 
 * 
 * @param t_f Highest boundary for time. CURRENTLY UNUSED AS FIXED TO MAX.
 * @param r_S Schwarzschild radius

 * @return map<double, double> : 
 * - integer key is lambdas' size -> value is number of such lambdas.
 * - double key n.5 is lambdas's size + 0.5 -> value is average \Delta r 
 * Also, note, we also keep:
 * - result[-1] = mintime;
 * - result[-2] = innermost;
 * - result[-3] = outermost;
 */
std::map<double,double> EmbeddedCauset::count_lambdas_withdr
(double& t_f, double r_S)
{
    if (strcmp(_spacetime._name, "BlackHole")==0)
    {
        if (_futures.size() == _size) /*if already defined*/
        {
            std::cout<<"No need for futs, straight to counting molecules\n";
            auto sizes = get_lambdas_sizes_withdr_from_futs(t_f,r_S);
            std::cout << "Finished get_lambdas_sizes" << std::endl;
            std::map<double,double> distr = get_lambdas_distr_withdr(sizes); 
            std::cout << "Finished get_lambdas_distr" << std::endl;
            return distr;
        }

        else if (_CMatrix.size()==0)
        {
            std::cout << "To create future, CMatrix must exist";
            throw std::invalid_argument("No CMatrix");
        }
        
        else
        {
            std::cout<<"Starting doing capped futs in count_lambdas"<<std::endl;
            _futures.resize(_size);
        
            #pragma omp parallel for
            for (int i=0; i<_size; i++)
            {
                int n_futs_of_i = 0;
                for (int j=i+1; j<_size; j++)
                {
                    #pragma omp critical
                    if (_CMatrix[i][j] == 1) {
                        _futures[i].insert(j);
                        n_futs_of_i += 1;
                    }

                    if (n_futs_of_i - 1 > 0){
                        break; /*break j loop, go to next i*/
                    }
                }
            }
            std::cout << "Finished done futs" << std::endl;
            auto sizes = get_lambdas_sizes_withdr_from_futs(t_f,r_S);
            std::cout << "Finished get_lambdas_sizes" << std::endl;
            std::map<double,double> distr = get_lambdas_distr_withdr(sizes); 
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
 * @brief   Finds lambdas in the causet connecting maximal elements 
 *          below t_f and inside the horizon with maximal-but-one elements
 *          outside the horizon.
 *          Currently works only for spacetime "Schwarzschild" in EForig coords, but
 *          could be expanded if needed in the future. 
 * 
 * @param t_f Highest boundary for time. CURRENTLY UNUSED AS FIXED TO MAX.
 * @param r_S Schwarzschild radius

 * @return map<double, pair<double,double>> : key is label of maximal element, 
 value.first is number of maximal but one elements associated with it,
 value.second is delta_r extension of molecule.
 */
std::map<int,std::pair<double, double> >
EmbeddedCauset::get_lambdas_sizes_withdr_from_futs (double& t_f, double r_S)
{
    if (!strcmp(_spacetime._name, "BlackHole")==0)
    {
        std::cout<<"Please choose 'BlackHole' for spacetime." <<
        "Other spacetimes might be available in the future" << std::endl;
        throw std::invalid_argument("Wrong spacetime");
    }

    // Maps label of maximal element 
    // to size and dr extension of its lambda
    std::map<int, std::pair<double, double>> lambdas;
    
    // To find point with lowest time component. Hypersurface set at t=tmax btw.
    std::vector<double> mintime_vec;
    std::vector<double> innermost_vec;
    std::vector<double> outermost_vec;
 
    for (int j = 1; j<_size; ++j)
    {
        // if j is maximal and inside the horizon
        if (_futures[j].size()==0 && _coords[j][1]<r_S) 
        {
            lambdas[j] = {0.,0.};
            for (int i = j-1; i>-1; --i)
            {
                if (_coords[j][0]>_coords[i][0]) //t_j>t_i SHOULD ALWAYS GO HERE
                {
                    // if i is maximal but one and outside the horizon
                    if (_futures[i].size()==1 && _coords[i][1]>r_S)
                    {
                        // if i-j is link
                        if (_futures[i].find(j) != _futures[i].end()) 
                        {
                            //update lambda count
                            lambdas[j].first += 1; 

                            // update dr if new i element further than previous
                            double dr_ij = _coords[i][1] - _coords[j][1];
                            if (dr_ij > lambdas[j].second)
                            lambdas[j].second = dr_ij * 1.;

                            // update statistics
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
    lambdas[-1] = {mintime, 0.};
    lambdas[-2] = {innermost, 0.};
    lambdas[-3] = {outermost, 0.};
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

 * @return std::map<double,double> : 
 * - int key is label of lambdas' size, value is number of occurrences.
 * - double key n.5 is lambdas' size + 0.5, value is average Delta r
 * - [-1] -> mintime
 * - [-2] -> innermost
 * - [.3] -> outermost
 */
std::map<double,double> EmbeddedCauset::get_lambdas_distr_withdr
        (const std::map<int, std::pair<double, double>> & lambdas)
{
    // Maps label of maximal element to size of its lambda
    std::map<double, double> lambdas_distr;

    for (auto pair : lambdas)
    {
        // pair.second = {size of j's molecule, dr of j's molecule}
        if (pair.first > 0 && pair.second.first != 0){
            lambdas_distr[pair.second.first] += 1;
            lambdas_distr[pair.second.first + 0.5] += pair.second.second;
        }

        else if (pair.first < 0)
            lambdas_distr[pair.first] = pair.second.first*1.;
    }

    for (auto pair : lambdas_distr)
    {
        // if pair is one of dr, hence pair.first = int.5
        // take average of drs, i.e. divide by number of size pair.firs mols
        if (pair.first > 0 && pair.first != (int)pair.first)
            lambdas_distr[pair.first] /= lambdas_distr[pair.first-0.5]; 
    }
    return lambdas_distr;
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