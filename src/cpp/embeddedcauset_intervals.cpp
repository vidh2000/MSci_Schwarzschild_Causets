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
 * @brief Get the interval between minimal and maximal element
 *          Essentially recreates the causet, getting rid of all
 *          the elements that aren't in the interval ->
 *          -> changes the pasts/futures sets and reduces the cmatrix
 * 
 *          The intervals are chosen at random and must contain #elements
 *          larger than min_size.
 * 
 * @param min_size : int. Minimal size of the interval (2 is default) 
 * @param max_size : int.  Max. size of the interval (size of causet is default)
 * @param N_max : int. Max number of tries to find the interval before stopping.
 *                (1000 is default).
 */
void EmbeddedCauset::get_random_intervals(int min_size, int max_size, int N_max) 
{

    if (min_size <=2){
        std::cout << "min_size>2 required!" << std::endl;
        throw std::runtime_error("");
    }

    bool found = false; 
    int N_tries = 0;

    if (max_size == 0)
    {
        max_size = _size;
    }
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
        if (n >= min_size && n<= max_size)
        {  

            // Interval includes a, b and the elements connecting them
            std::unordered_set<int> interval = set_intersection(
                        _futures[a], _pasts[b]);
            interval.insert(a);
            interval.insert(b);

            // Find indices to remove i.e all but the inclusive interval
            std::unordered_set<int> indices_to_remove = set_diff(
                                    all_indices,interval);
                    
            // Create a sorted vector of indices to discard
            std::vector<int> to_discard(indices_to_remove.begin(),
                                        indices_to_remove.end());
            std::sort(to_discard.begin(),to_discard.end());
            
            // Create a sorted vector of remaining indices (the interval)
            std::vector<int> ordered_interval(interval.begin(),
                                              interval.end());
            std::sort(ordered_interval.begin(),ordered_interval.end());              
            

            // Assumes matrix is created and future and pasts but no links.
            EmbeddedCauset::discard(to_discard,ordered_interval,
                                                true,true,false);

            found = true;
        }
        else{
            N_tries +=1;
            continue;
        }
    }
}


/**
 * @brief Get the average number of chains of size k, up to 
 *        including size k_max,
 *        in an interval of size (min_size and max_size) in a causet
 *        over "N_intervals" random intervals.
 *        REQUIRES CMATRIX, AND PAST AND FUTURE SETS 
 * 
 * @param N_intervals - Number of intervals
 * @param min_size - Minimal size of the interval (min # of elements in it) 
 * @param k_max - maximal (included) length of chain we care about
 * @param max_size - Max. size of the interval (max # of elements in it
 *                                              == _size by default)
 * @param N_max - max number of tries to find the interval before stopping
 * @param avoid_boundaries bool : True (Not-default) implies that the extremi
 * of the intervals are within 25 and 75 % of space interval to avoid boundary
 * effects.
 * It assumes coords[1] is a radial distance from origin.
 * 
 * @exception std::runtime_error - if does not find suitable interval in N_max 
 * tries.
 * 
 * @return Returns a vector of length N_intervals:
 *         - for each interval there's a pair <N_chains_k,r_avg>
 *              - N_chains_k: 
 *                  a vector with "k_max" entries for chains of length 1..k    
 *                  storing the number of such chains in the interval
 *              - r_avg:     
 *                  an average r value of the interval
 *              
 *         i.e for each interval, for each chain of length k
 *          you have information about the number of such chains
 *          and what was the average r-value for the interval
 *          (to be able to check if it varies w.r.t r in Schwarzschild) 
 */
vector<std::pair<vector<double>,double>> EmbeddedCauset::get_Nchains_inInterval(
                    int N_intervals, int min_size, int k_max,
                    int max_size, int N_max, bool avoid_boundaries)
                     //==0, 1000
{
    if (min_size <=2){
        std::cout << "min_size>2 required!" << std::endl;
        throw std::runtime_error("");
    }

    if (max_size == 0)
    {
        max_size = _size;
    }

    // Define vars and outcome vars
    int N_intervals_found = 0;
    vector<std::pair<vector<double>,double>> results;

    // Define limits if avoid_boundaries
    double rmin, rmax;
    if (avoid_boundaries)
    {
        std::vector<double> center = _shape._center; 
        double duration = _shape._params.find("duration")->second;
        double radius   = _shape._params.find( "radius" )->second;
        double hollow   = _shape._params.find( "hollow" )->second;
        rmin = (hollow != 0.)? radius*hollow + 0.25*radius*(1-hollow) : 0.;
        rmax = 0.75*radius;
    }
    

    while (N_intervals_found<N_intervals)
    {
        int N_tries = 0;
        bool found = false; 
        vector<double> chain_arr;
        double r_avg = 0;

        while (!found)
        {
            // Failsafe
            if (N_tries > N_max){
                std::cout << "Couldn't find suitable interval in " << N_max
                    << "tries" << std::endl;
                throw std::runtime_error("");
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
            else if (e1<e2){
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

            // to have an interval require a prec b, so skip if not
            if (_CMatrix[a][b] == 0){
                N_tries += 1;
                continue;
            }

            // Check they respect boundaries, if that ise demanded
            if (avoid_boundaries){
                double r1 = _coords[e1][1];
                double r2 = _coords[e2][1];
                if (!(rmin <= r1 && r1 <= rmax && rmin <= r2 && r2 <= rmax)){
                    N_tries += 1;
                    continue;
                }
            }

            int n = IntervalCard(a, b);
            if (n >= min_size && n<= max_size)
            {  

                /*          //OLD WAY WITH SET INTERSECTION

                // Create set_intersection for cmatrix..
                // Interval includes a, b and the elements connecting them
                std::unordered_set<int> interval = set_intersection(
                            _futures[a], _pasts[b]);
                interval.insert(a);
                interval.insert(b);

                // Create a sorted vector of the interval (remaining indices)
                std::vector<int> ordered_interval(interval.begin(),
                                                interval.end());
                std::sort(ordered_interval.begin(),ordered_interval.end());              
                */
                

                // Create set_intersection for cmatrix..
                // Interval includes a, b and the elements connecting them
                vector<int> ordered_interval = {a};
                for (int i=a+1; i<b; i++) {
                    if(_CMatrix[a][i] && _CMatrix[i][b]) {
                        ordered_interval.push_back(i);
                    }
                }
                ordered_interval.push_back(b);

                if(ordered_interval.size() != n) {
                    std::cout << "n="<<n<<", interval size="<<
                    ordered_interval.size()<<std::endl;
                }

                // Create a copy of the (cut) interval "reduced" cmatrix
                vector<vector<int>> M =
                            EmbeddedCauset::getIntervalCmatrix(ordered_interval);
                
                // Get the number of chains up to including size k
                // where k=1 == _size
                double C1 = (double)n;
                double C2 = sumMatrix(M);

               if (k_max >=3) { 
                    vector<vector<int>> M2 = matmul(M,M);
                    double C3 = sumMatrix(M2);
                    chain_arr.push_back(C1);
                    chain_arr.push_back(C2);
                    chain_arr.push_back(C3);

                    if (k_max == 4){
                        vector<vector<int>> M3 = matmul(M2,M);
                        double C4 = sumMatrix(M3);
                        chain_arr.push_back(C4);
                    }
               }
                else if (k_max>4) {
                    print("Haven't implemented this for k>4!");
                    throw std::runtime_error("");
                }
                else {
                    print("What did you choose for k?");
                    print("k<=4 required!");
                    throw std::runtime_error("");
                }

                // Find the average r value in the interval
                for (int index : ordered_interval) {
                    r_avg += _coords[index][1];
                }
                r_avg = r_avg/(double)n;
        

                // Found the suitable interval in the causet
                found = true;

                //Create the pair of <N_chainK vector, r_avg>
                std::pair<vector<double>,double> interval_result =
                                                        {chain_arr,r_avg};
                results.push_back(interval_result);

                N_intervals_found++;
                std::cout << "Found "<<
                    (N_intervals_found) << "/" << N_intervals <<
                    " intervals" << std::endl;
            }
            else{
                N_tries +=1;
                continue;
            }
        }   
    }
    // Found number of (1...k_max)-sized chains for N_intervals
    return results;  
}