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
#include <random>

#include "functions.h"
#include "MyVecFunctions.h"
#include "causet.h"

using std::vector;
using std::set;
using std::unordered_set;


/**
 * @brief Causet class.
 *
 * 
 * @param causet: a vector of vectors of integers. //not implemented
 *  Essentially it is the causal matrix where (M)_{i,j}
 *  can take value of 1 if e_j<e_i, 0 if they aren't related and
 *  -1 if e_j<e_i and they are a link.
 * @param pasts: a vector of sets for all elements.
 * 
 * @param pastLinks: a vector of sets containing past links to each element
 */

// CONSTRUCTORS
Causet::Causet(){}
Causet::Causet(vector<vector<double>> Cmatrix, 
                bool past_links, // = false,
                bool fut_links) // = false);
{

}
Causet::Causet(vector<vector<double>> coordinates,
               const char* method) // = "pasts");
{

}

Causet::Causet(vector<vector<double>> coordinates,
                const char* method)
{

}

///////////////////////////////////////////////////////////////////////////////
// ORDERING FRACTION FUNCTIONS (OVERRIDING FOR TYPES)
///////////////////////////////////////////////////////////////////////////////
/**
 * @brief   Find the ordering fraction of an interval: 
            The ratio of actual relations over possible ones in such interval.
 * 
 * @param mode: string
            Use as denominator:
            - 'choose' -> |A|(|A|-1)/2, i.e. |A| choose 2 (Default).
            - 'n2'     -> (|A|^2)/2.
 *            
 * @return  Ordering fraction of Alexandrov Interval
            This is nrelations / (N choose 2)
 */

/**
 * @brief   Ordering fraction determined for a "Causet" object.
 *          Can either be determined via a matrix or via pasts/futures sets.
 *          See above description above all ord_fr definitions.
 */
double Causet::ord_fr(Causet A,
            const char* denominator, // = "choose",
            bool isdisjoined) // = true);
{
    if (A._CMatrix.size())
    {
        return Causet::ord_fr(A._CMatrix, denominator, isdisjoined);    
    }
    else if (A._pasts.size())
    {
        return Causet::ord_fr(A._futures,A._pasts,denominator,isdisjoined);
    }
}

/**
 * @brief   Ordering fraction determined from a matrix.
 *          See above description above all ord_fr definitions
 */
double Causet::ord_fr(vector<vector<int8_t>> A,
                const char* denominator,// = "choose",
                bool isdisjoined)// = true);
{

}

/**
 * @brief   Ordering fraction determined from sets of futures and pasts.
 *          See above description above all ord_fr definitions
 */
template<typename SET>
double Causet::ord_fr(vector<SET> A_futures, 
                vector<SET> A_pasts,
                const char* denominator,// = "choose",
                bool isdisjoined)// = true);
{
    if (_CMatrix.size())
    {

    }
    else if (_pasts.size())
    {
        
    }
}

/**
 * @brief   Find the ordering fraction of an interval between a and b: 
            the ratio of actual relations over possible in such interval.
 * 
 * @param a: integer "causetevent", a<b wanted
 * @param b: integer "causetevent", a<b wanted
 */
double Causet::ord_fr(int a, int b,
                const char* denominator)// = "choose"
{

}

///////////////////////////////////////////////////////////////////////////////
// Dimension estimator
///////////////////////////////////////////////////////////////////////////////
/**
 * @brief Use Myrheim-Meyers dimensional estimator to compute the 
          fractal dimension (not necesseraly int).
          Only intended to work with a Causet with
          _pasts, and _futures vector<set<int>> objects existing.
        
 * 
 * @param method: str 
            - 'random': randomly sample.
            - 'big': take all events with no past, all with no future
                     and apply estimator ro their combinations.
 *  
 * @param d0: float 
            Initial guess for dimension.
            Default is 2.
 * 
 * @param Nsamples: int 
            Times to iterate procedure to then average on if method "random".
            Default is 20.
 * 
 * @param size_min: int\n
            Minimum size of Alexandrov Sets on which you apply estimators.
            Default is 20 elements.
 * 
 * @param size_max: int\n
            Maximum size of Alexandrov Sets on which you apply estimators.
            Default and highly recommended is np.Inf.
 * 
 * @return
        - dimension estimate: float
        - dimension std: float
 */
double Causet::MM_drelation(double d)
{
    double a = std::tgamma(d+1);
    double b = std::tgamma(d/2);
    double c = 4* std::tgamma(3*d/2);
    return a*b/c;
}


vector<double> Causet::MMdim_est(const char* method,// = "random",
                                int d0,// = 2,
                                int Nsamples,// = 20,
                                int size_min,// = 10,
                                double size_max)// = nan("")
{
    std::cout << "NOTE: MMd works only in flat spacetime" << std::endl;

    auto MM_to_solve = [](double d, double ord_fr){
        return Causet::MM_drelation(d) - ord_fr/2;
        };
    
    // Variables to be used
    int* N = &_size;
    std::vector<double> destimates;
    
    if (method == "random")
    {
        int isample = 0;
        int fails = 0;
        int successes = 0;

        while (isample < Nsamples)
        {
            if ((fails>= 1000) && (successes == 0))
            {
                std::cout << "Found 0/1000 OK Alexandrov intervals. \
                Causet portion too smol. Returning Dim<0 values.";
                return {-1,-1};
            }

            // Pick two random elements

            // Define mersenne_twister_engine Random Gen. (with random seed)
            std::random_device rd;
            int seed = rd();
            std::mt19937 gen(seed);
            std::uniform_real_distribution<> dis(0,*N);
            int e1 = (int) dis(gen), e2 =(int) dis(gen);
            int a;
            int b;
            if (e1 == e2){
                fails += 1;
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
                fails += 1;
                continue;
            }
            
            int n = IntervalCard(a, b);
            if (n >= size_min && n<= size_max) 
            {
                successes += 1;
                double fr_i = this->ord_fr(a,b,"choose");
                if (fr_i ==1)
                {
                    destimates.push_back(1);
                    isample +=1;
                }
                else
                {
                    double fr_i = fr_i * (n-1)/n; //correction for MMestimator
                    double min = 0;
                    double max = 5;
                    // Estimate dimension of Causet
                    double d_i = bisect(min,max,......);
                    destimates.push_back(d_i);
                    isample +=1;
                }
            }
            else
            {
                fails +=1;
                continue;
            }
        }
    }
    else if (method == "big")
    {
        vector<int> As = {};
        vector<int> Bs = {};
        for (int e = 0; e<*N; e++)
        {
            if (_pasts[e].size() == 0){
                As.push_back(e);
            }
            else if (_futures[e].size() == 0){
                Bs.push_back(e);
            }
        }
        for (int i=0; i<As.size(); i++)
        {
            for (int j=0; j<Bs.size(); j++)
            {   
                int a = As[i];
                int b = Bs[j];
                int n = IntervalCard(a, b);
                if (n >= size_min && n<= size_max) 
                {
                    double fr_i = this->ord_fr(a,b,"choose");
                    if (fr_i ==1)
                    {
                        destimates.push_back(1);
                    }
                    else
                    {
                        double fr_i = fr_i * (n-1)/n; //correction for MMestimator
                        double min = 0;
                        double max = 5;
                        // Estimate dimension of Causet
                        double d_i = bisect(min,max,......);
                        destimates.push_back(d_i);
                    }
                }
            }
        }
    }
    else
    {
        const char* errormsg = "Choose 'method' parameter 'random' or 'big'";
        throw std::invalid_argument(errormsg);
    }

    // Return mean and std of dimension estimate result
    double sum = std::accumulate(std::begin(destimates),
                                 std::end(destimates), 0.0);
    double mean = sum/destimates.size();

    double accum = 0.0;
    std::for_each(std::begin(destimates), std::end(destimates),
                    [&](const double d)
                    {
                    accum += (d - mean) * (d - mean);
                    });

    double stdev = sqrt(accum / (destimates.size()));
    vector<double> result = {mean,stdev};
    return result;
    
}   





//=============================================================================
//=============================================================================
//MODIFIERS   //===============================================================
//=============================================================================
//=============================================================================
/**
 * @brief Coarse grain Causet of "card" events.
 * 
 * @param card : int
 * @param make_matrix : if true and _CMatrix non-empty, update _CMatrix 
 * @param make_sets : if true, update _pasts and/or _futures, if they are 
 * defined
 * @param make_links : if true, update _past and/or _future links, if they are 
 * defined
 */
void Causet::coarsegrain(int card, bool make_matrix, 
                         bool make_sets, bool make_links)
{
    vector<int> labels = Causet::distinct_int_random(card, _size);
    this->discard(labels, make_matrix, make_sets, make_links); 
}
void Causet::cgrain(int card, bool make_matrix, 
                    bool make_sets, bool make_links)
{
    vector<int> labels = Causet::distinct_int_random(card, _size);
    this->discard(labels, make_matrix, make_sets, make_links); 
}
void Causet::coarsegrain(double fract, bool make_matrix, 
                         bool make_sets, bool make_links)
{
    int card = fract * _size;
    vector<int> labels = Causet::distinct_int_random(card, _size);
    this->discard(labels, make_matrix, make_sets, make_links); 
}
void Causet::cgrain(double fract, bool make_matrix, 
                    bool make_sets, bool make_links)
{
    int card = fract * _size;
    vector<int> labels = Causet::distinct_int_random(card, _size);
    this->discard(labels, make_matrix, make_sets, make_links); 
}
/**
 * @brief Discard event number "label" (from 0 to N-1)
 * 
 * @param label : int
 * @param make_matrix : if true, update _CMatrix, if already defned
 * @param make_sets : if true, update _pasts and/or _futures, if they are 
 * defined
 * @param make_links : if true, update _past and/or _future links, if they are 
 * defined
 */
void Causet::discard(int label, bool make_matrix = true, 
                bool make_sets = false, bool make_links = true)
{
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


/**
 * @brief Discard events numbered by i in "label" (from 0 to N-1)
 * 
 * @param label : vector<int>
 * @param make_matrix : if true, update _CMatrix, if already defned
 * @param make_sets : if true, update _pasts and/or _futures, if they are 
 * defined
 * @param make_links : if true, update _past and/or _future links, if they are 
 * defined
 */
void Causet::discard(vector<int> labels, bool make_matrix = true, 
                bool make_sets = false, bool make_links = true)
{
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


void Causet::discard_from_set(unordered_set<int> &myset, int label)
{
    int N = myset.size();
    unordered_set<int> buffer;
    for (int j : myset)
    {
        if (j<label)
            {buffer.insert(j);}
        if (j>label)
            {buffer.insert(j-1);}
    }
    myset = buffer;
}

template<typename m>
void Causet::discard_from_set(unordered_set<m> &myset, vector<int> labels)
{
    labels.sort();
    int N = myset.size();
    unordered_set<m> buffer;
    int startpoint = 0;
    for (m j : myset) //not ordered
    {
        if (j>labels[-1])
            {buffer.insert[j-labels.size()];}
        else
        {
            for (int s = 0; s<labels.size(); s++)
                {
                    if (labels[s] == j)
                        {break;}
                    else if (labels[s] > j)
                        {buffer.insert[j-s]; break;}
                }
        }
    }
    myset = buffer;
}


/**
 * @brief Generate "size" random ints between 0 and N-1.
 */
vector<int> Causet::distinct_int_random(int size, int N, int seed)
{
    if (size < N/2) {return distinct_int_random1(size, N, seed);}
    else {return distinct_int_random2(size, N, seed);}
}


/**
 * @brief Generate "size" random ints between 0 and N-1.
 */
vector<int> Causet::distinct_int_random1(int size, int N, int seed)
{
    vector<int> result(size);
    if (!seed)
        {
         auto rd = std::random_device {};
         seed = rd();
        }
    for(int i = 0; i < size; ++i)
    {
    int r;
    while(std::find(result.begin(), result.end(), r) != result.end())
    {
        r = rand(seed)%N;
        result[i] = r;   
    }
    return result;
}

/**
 * @brief Generate "size" random ints between 0 and N-1.
 */
vector<int> Causet::distinct_int_random2(int size, int N, int seed)
{
    vector<int> result(N);
    if (!seed)
        {
         auto rd = std::random_device {};
         seed = rd();
        }
    for(int i = 0; i < N; ++i)
        {result[i] = i;}
    auto rng = std::default_random_engine{seed};

    std::shuffle(std::begin(resul), std::end(result), rng);
    result = result.resize(size);
    return result;
}
