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
#include "vecfunctions.h"
#include "causet.h"

using std::vector;
using std::set;
using std::unordered_set;


/**
 * @brief Causet class.
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
{}

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
            const char* denominator) // = "choose"
{
    if (A._CMatrix.size())
    {
        return Causet::ord_fr(A._CMatrix, denominator);    
    }
    else if (A._pasts.size())
    {
        return Causet::ord_fr(A._pasts,denominator);
    }
    else{
        const char* errormsg = "Provided Causet is empty! Size matters...";
        std::cout << errormsg << " Returning ord_fr = 0.0" << std::endl;
        return 0.0;
    }
}

/**
 * @brief   Ordering fraction determined from a CMatrix.
 *          See above description above all ord_fr definitions
 */
double Causet::ord_fr(vector<vector<int8_t>> M,
                const char* denominator)// = "choose",
{
    if (denominator!= "choose" || denominator!= "n2"){
        throw std::invalid_argument("Param 'denominator' must be \
                                'choose' or 'n2'");}
    int N = M.size();
    int nrelations = 0;
    for (int j; j<N; j++){
        for (int i; i<j;j++){
            if (M[i][j] != 0){
                nrelations +=1;
            }
        }
    }
    double fr = 2 * nrelations/ (N* (N - (denominator=="choose")));
    return fr;
}

/**
 * @brief   Ordering fraction determined from sets of futures and pasts.
 *          See above description above all ord_fr definitions
 */
template<typename SET>
double Causet::ord_fr(vector<SET> A_pasts,
                const char* denominator)// = "choose",
{
    if (denominator!= "choose" || denominator!= "n2"){
        throw std::invalid_argument("Param 'denominator' must be \
                                'choose' or 'n2'");}
    int N = A_pasts.size();
    int nrelations = 0;
    for (auto e_i : A_pasts)
    {
        nrelations += e_i.size();
    }
    double fr = 2 * nrelations/ (N* (N - (denominator=="choose")));
    return fr;
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
    if (denominator!= "choose" || denominator!= "n2"){
        throw std::invalid_argument("Param 'denominator' must be \
                                'choose' or 'n2'");}
     
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
    /*
    Dimension function to be solved for with
    order fraction for MM-dimension estimation.
    */
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

    
    // Variables to be used
    int* N = &_size;
    vector<double> destimates;
    

    if (method == "random")
    {
        int fails = 0;
        int successes = 0;
        while (Nsamples>0)
        {
            if (fails>= 1000 && successes == 0)
            {
                std::cout << "Found 0/1000 OK Alexandrov intervals. \
                Causet portion too small. Returning {-1,-1} values.";
                vector<double> returnerr = {-1,-1};
                return returnerr;
            }

            // Pick two random elements
            

            // Define mersenne_twister_engine Random Gen. (with random seed)
            std::random_device rd;
            int seed = rd();
            std::mt19937 gen(seed);
            std::uniform_real_distribution<> dis(0,*N);
            int e1 = (int) dis(gen), e2 =(int) dis(gen);
            int a; int b;
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
                    Nsamples --;
                }
                else
                {
                    //Order fraction correction for MMestimator
                    double fr_i = fr_i * (n-1)/n;

                    // Define function whose root needs to be found
                    auto MM_to_solve = [fr_i](double d){
                        return Causet::MM_drelation(d) - fr_i/2;};

                    double dmin = 0.75;
                    double dmax = 10;
                    // Estimate dimension of Causet
                    double d_i = bisection(MM_to_solve,dmin,dmax);
                    destimates.push_back(d_i);
                    Nsamples --;
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
                        //Order fraction correction for MMestimator
                        double fr_i = fr_i * (n-1)/n;

                        // Define function whose root needs to be found
                        auto MM_to_solve = [fr_i](double d){
                            return Causet::MM_drelation(d) - fr_i/2;};

                        double dmin = 0.1;
                        double dmax = 5;
                        // Estimate dimension of Causet
                        double d_i = bisection(MM_to_solve,dmin,dmax);
                        destimates.push_back(d_i);
                        Nsamples --;
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
//REP & SAVE  //===============================================================
//=============================================================================
//=============================================================================
void Causet::saveC(const char* path_file_ext)
{
    std::ofstream out;
    out.open(path_file_ext);
    for (auto row : _CMatrix) 
    {
        for (auto col : row)
            {out << col <<',';}
        out<<std::endl;
    }
    out.close();
    return;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//\\ INTERVALS   \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/**
 * @brief Compute cardinality of causality interval between a and b.
 * 
 * @param a int : label of event a.
 * @param b int : label of event b
 * @param includeBoundary bool : innclude a and b in count?
 * @return int : Cardinality of Interval
 */
int Causet::IntervalCard(int a, int b, bool includeBoundary)
{
    if (a==b)
        {return 1;}
    // IF DEFINED USE CMATRIX
    if (_CMatrix.size())
    {
        int Nintersections = 2 * includeBoundary;
        if (a<b)
        {
            for (int i = a+1; i<b; i++)
                {Nintersections += _CMatrix[a][i] && _CMatrix[i][b];}
        }
        else
        {
            for (int i = b+1; i<a; i++)
                {Nintersections += _CMatrix[i][a] && _CMatrix[b][i];}
        }
        return Nintersections;
    }
    // ELSE USE SETS
    else //if ( _pasts.size() && _futures.size())
    {
        if (a<b)
        {
            int Nintersections = 2 * includeBoundary;
            if (_futures[a].size()<_pasts[b].size()) //loop over shortest
            {
                for (int e_ai : _futures[a])
                    {Nintersections += _pasts[b].find(e_ai) !=_pasts[b].end();}
            }
            else
            {
                for (int e_bi : _pasts[b])
                    {Nintersections+= _futures[a].find(e_bi)!=_futures[a].end();}
            }
            return Nintersections;
        }
        else /*b<a*/
        {
            int Nintersections = 2 * includeBoundary;
            if (_pasts[a].size()<_futures[b].size()) //loop over shortest
            {
                for (int e_ai : _pasts[a])
                    {Nintersections+= _futures[b].find(e_ai)!=_futures[b].end();}
            }
            else
            {
                for (int e_bi : _futures[b])
                    {Nintersections += _pasts[a].find(e_bi) !=_pasts[a].end();}
            }
            return Nintersections;
        }
    }  
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
                         bool make_sets, bool make_links, int seed)
{
    vector<int> labels = distinct_randint(card, _size, seed);
    this->discard(labels, make_matrix, make_sets, make_links); 
}
void Causet::cgrain(int card, bool make_matrix, 
                    bool make_sets, bool make_links, int seed)
{
    vector<int> labels = distinct_randint(card, _size, seed);
    this->discard(labels, make_matrix, make_sets, make_links); 
}

/**
 * @brief Coarse grain Causet of "fract*_size" events.
 * 
 * @param fract : double, fraction of elements to remove
 * @param make_matrix bool: if true and _CMatrix non-empty, update _CMatrix 
 * @param make_sets bool: if true, update _pasts and/or _futures, if they are 
 * defined
 * @param make_links bool: if true, update _past and/or _future links, if they 
 * are defined
 */
void Causet::coarsegrain(double fract, bool make_matrix, 
                         bool make_sets, bool make_links, int seed)
{
    int card = fract * _size;
    vector<int> labels = distinct_randint(card, _size, seed);
    this->discard(labels, make_matrix, make_sets, make_links); 
}
void Causet::cgrain(double fract, bool make_matrix, 
                    bool make_sets, bool make_links, int seed)
{
    int card = fract * _size;
    vector<int> labels = distinct_randint(card, _size, seed);
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
void Causet::discard(int label, bool make_matrix, 
                     bool make_sets, bool make_links)
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
void Causet::discard(vector<int> labels, bool make_matrix, 
                    bool make_sets, bool make_links)
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

Causet::~Causet() {}

int main(){

    std::cout << "causet.cpp WORKS!" << std::endl;
}

