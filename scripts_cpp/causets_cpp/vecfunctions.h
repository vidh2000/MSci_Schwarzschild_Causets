#ifndef MYCFUNCTIONS_H
#define MYCFUNCTIONS_H


#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <random>
#include <set>
#include <stack>
#include <stdio.h>
#include <stdlib.h>
#include <stdexcept>
#include <string>
#include <vector>
#include <omp.h>
#include "functions.h"

//_____________________________________________________________________
//
//---------------------------------------------------------------------
// RANDOM GENERATORS FUNCTIONS
//---------------------------------------------------------------------
//_____________________________________________________________________
/*
 *  Generates a matrix of size (rows, cols) with
 *  random entries lying in range [min,max]
*/
inline
std::vector<std::vector<double>> generate_2Dvector(int rows, int cols, 
                                         double min, double max)
{
    srand(time(0)); // to generate pseudo random numbers each time
    std::vector<std::vector<double>> matrix;
    matrix.resize(rows);
    
    // Random generator stuff
    std::random_device rd;  
    std::mt19937 gen(rd()); 
    std::uniform_real_distribution<> dis(min,max);

    for(int i=0; i<rows; i++)
    {
        matrix[i].resize(cols);
        for (int j=0; j<cols; j++){
            matrix[i][j] = dis(gen);
        }
    }
    return matrix;
}

inline
std::vector<std::vector<int>> generate_2Dvector_of_ints(int rows, int cols, 
                                         int min, int max)
{

    
    std::vector<std::vector<int>> matrix;
    matrix.resize(rows);
    
    std::srand(time(0)); //Randomise seed initialisation
    
    for(int i=0; i<rows; i++)
    {
        matrix[i].resize(cols);
        for (int j=0; j<cols; j++){
            matrix[i][j] = rand() % 2;
        }
    }
    return matrix;
}


/**
 * @brief Generate "size" DIFFERENT random ints between 0 and N-1.
 */
inline
std::vector<int> distinct_randint1(int size, int N, int seed=0)
{
    std::vector<int> result;
    result.resize(size);

    if (seed==0){
        std::random_device rd;
        seed = rd();}
    std::cout<<seed<<std::endl;

    std::mt19937 gen(seed);
    std::uniform_real_distribution<> dis(0, N);
    
    for(int i = 0; i < size; ++i)
    {
        int r;
        while(std::find(result.begin(), result.end(), r) != result.end())
        {
            r = (int) dis(gen);  
        }
        result[i] = r; 
    }
    return result;
}

/**
 * @brief Generate "size" DIFFERENT random ints between 0 and N-1.
 */
inline
std::vector<int> distinct_randint2(int size, int N, int seed=0)
{
    std::vector<int> result(N); 
    if (seed==0)
    {
         std::random_device rd;
         seed = rd();
    }
    for(int i = 0; i < N; ++i)
        {result[i] = i;}
    
    // Narrow seed int to fit into range [0,2147483647]
    seed = abs(seed);
    auto rng = std::default_random_engine(seed);

    std::shuffle(std::begin(result), std::end(result), rng);
    result.resize(size);
    return result;
}

/**
 * @brief Generate "size" DIFFERENT random ints between 0 and N-1.
 */
inline
std::vector<int> distinct_randint(int size, int N, int seed = 0)
{
    //RANDINT1 NOT WORKING AT THE MOMENT
    if (size < 0 ){
        return distinct_randint1(size, N, seed);}
    else{
        return distinct_randint2(size, N, seed);}
}



//_____________________________________________________________________
//
//---------------------------------------------------------------------
// CLASSIC FUNCTIONS
//---------------------------------------------------------------------
//_____________________________________________________________________
//From https://thispointer.com/cpp-vector-print-all-elements/

/** @brief Prints given vector, with given separator.
*
*  @param vec The vector to print. Type: std::vector<whatever>.
*  @param sep Separator to use between entries. Type: Char. Default " , ".
*  @param brackets Start and end with brackets?. Type: bool. Default true.
* 
*  @return Void.
**/
template <typename T>
inline
void print_vector(const std::vector<T> & vec, 
                  std::string sep=", ",
                  bool brackets = true)
    {   
        int size = vec.size();
        if (brackets) {std::cout<<"{";}
        for(int index = 0; index < size ; index++)
        {   
            if (index == 0) {std::cout<<vec[index];}
            else {std::cout<<sep<<vec[index];}
        }
        if (brackets) {std::cout<<"}";}
        std::cout<<std::endl;  
}

/** @brief Prints given vector, with given separator.
*
*  @param vec The vector to print. Type: std::vector<whatever>.
*  @param sep Separator to use between entries. Type: Char. Default " , ".
*  @param brackets Start and end with brackets?. Type: bool. Default true.
* 
*  @return Void.
**/
template <typename T>
inline
void print(const std::vector<T> & vec, 
                  std::string sep=", ",
                  bool brackets = true)
    {   
        int size = vec.size();
        if (brackets) {std::cout<<"{";}
        for(int index = 0; index < size ; index++)
        {   
            if (index == 0) {std::cout<<vec[index];}
            else {std::cout<<sep<<vec[index];}
        }
        if (brackets) {std::cout<<"}";}
        std::cout<<std::endl;  
}


template <typename T>
inline
void print_vector(std::vector<std::vector<T>> vec, 
                  std::string sep=", ",
                  bool brackets = true)
    /** @brief rints given vector, with given separator.
    *
    *  @param vec The vector to print. Type: std::vector<vector<whatever>>.
    *  @param sep Separator to use between entries. Type: Char. Default " , ".
    *  @param brackets Start and end with brackets?. Type: bool. Default true.
    * 
    *  @return Void.
    **/
    {   int size = vec.size();
        if (brackets) {std::cout<<"{";}
        for(int index = 0; index < size ; index++)
        {
            std::vector<T> vind = vec[index];
            if (index != 0)
            {std::cout<<" {";}
            else {std::cout << "{";}
            for (int j = 0; j < vind.size(); j++)
            {
                if (j == 0) {std::cout<<vind[0];}
                else {std::cout<<sep<<vind[j];}
            }
            std::cout << "}";

            if (index != size -1) {std::cout<< sep << std::endl;}
        }
        if (brackets) {std::cout<<"}";}
        std::cout<<std::endl;  
}


template <typename T>
inline
void print(std::vector<std::vector<T>> vec, 
                  std::string sep=", ",
                  bool brackets = true)
    /** @brief rints given vector, with given separator.
    *
    *  @param vec The vector to print. Type: std::vector<vector<whatever>>.
    *  @param sep Separator to use between entries. Type: Char. Default " , ".
    *  @param brackets Start and end with brackets?. Type: bool. Default true.
    * 
    *  @return Void.
    **/
    {   int size = vec.size();
        if (brackets) {std::cout<<"{";}
        for(int index = 0; index < size ; index++)
        {
            std::vector<T> vind = vec[index];
            if (index != 0)
            {std::cout<<" {";}
            else {std::cout << "{";}
            for (int j = 0; j < vind.size(); j++)
            {
                if (j == 0) {std::cout<<vind[0];}
                else {std::cout<<sep<<vind[j];}
            }
            std::cout << "}";

            if (index != size -1) {std::cout<< sep << std::endl;}
        }
        if (brackets) {std::cout<<"}";}
        std::cout<<std::endl;  
}



template <typename T>
inline
T myvecsum(std::vector <T> v)
    {return std::accumulate(v.begin(),v.end(), .0);}

template <typename T, typename F>
inline
double myvecsum(std::vector <T> v, F func)
{
    double sum = 0;
    for (int i = 0; i < v.size(); i++) sum += func(v[i]);
    return sum;
}


template <typename T1>
inline
T1 vecmax(std::vector<T1> v)
    {return *(std::max_element(v.begin(), v.end()));}

template <typename T1>
inline
T1 vecmin(std::vector<T1> v)
    {return *(std::min_element(v.begin(), v.end()));}


template <typename T1, typename T2>
int argmax(std::vector<T1, T2> const& v, int begin = 0, int end = 0) 
/** @brief Gets index of maximum in vector in interval [begin, end)
 *
 *  @param v     Vector where to look in.
 *  @param begin Starting point of range. Type int. 
 *               Default 0.
 *  @param end   Endng point of range. Type int.
 *               Default 0 (meaning v.end() is end of range). 
 * 
 *  @return Index of maximum.
 **/
{
    auto b = v.begin() + begin;
    auto e = (end == 0) ? v.end() : v.begin()+end;
    return static_cast<int>(std::distance(v.begin(), 
                                          max_element(b, e)));
}

template <typename T1, typename T2>
int argmin(std::vector<T1, T2> const& v, int begin = 0, int end = 0) 
{
/** @brief Gets index of minimum in vector in interval [begin, end)
 *
 *  @param v     Vector where to look in.
 *  @param begin Starting point of range. Type int. 
 *               Default 0.
 *  @param end   Endng point of range. Type int.
 *               Default 0 (meaning v.end() is end of range). 
 * 
 *  @return Index of minimum.
 **/
    auto b = v.begin() + begin;
    auto e = (end == 0) ? v.end() : v.begin()+end;
    return static_cast<int>(std::distance(v.begin(), 
                                          min_element(b, e)));
}


/**
 * @brief Produces a copy of a matrix which was cut to
 *        only rows and columns that correspond to the indices -
 *        elements from the interval
 * @param matrix :  NxN matrix
 * @param interval: vector which contains the target indices telling which
 *                  rows and columns not to cut away
 */
template <typename T>
inline
std::vector<std::vector<T>> get_reducedMatrix(std::vector<std::vector<T>> &matrix,
                std::vector<int> interval)
{

    // Make sure values in the interval are ordered
    std::sort(interval.begin(),interval.end());
    
    int N = matrix.size();
    std::vector<std::vector<T>> new_matrix;

    for (int i = 0; i<N; i++)
    {
        if (std::find(interval.begin(),interval.end(),i) == interval.end()) {
            continue;
        }
        std::vector<T> row;
        for (int j=0; j<N; j++)
        {
            if (std::find(interval.begin(),interval.end(),j) == interval.end()) {
            continue;
            }

            row.push_back(matrix[i][j]);        
        }
        new_matrix.push_back(row);
    }
    return new_matrix;
}



/**
 * @brief Cuts the Cmatrix to only rows and columns that correspond to
 *          the indices - elements from the interval
 * @param matrix : i.e CMatrix 
 * @param interval: vector which contains the target indices telling which
 *                  rows and columns not to cut away
 */
template <typename T>
inline
void remove_indices_fromCmatrix(std::vector<std::vector<T>> &matrix,
                std::vector<int> interval)
{

    // Make sure values in the interval are ordered
    std::sort(interval.begin(),interval.end());
    
    int N = matrix.size();
    std::vector<std::vector<T>> new_matrix;

    for (int i = 0; i<N; i++)
    {
        if (std::find(interval.begin(),interval.end(),i) == interval.end()) {
            continue;
        }
        std::vector<T> row;
        for (int j=0; j<N; j++)
        {
            if (std::find(interval.begin(),interval.end(),j) == interval.end()) {
            continue;
            }

            row.push_back(matrix[i][j]);        
        }
        new_matrix.push_back(row);
    }
    matrix = new_matrix;
}


////// idk that doesn't really work, still keeping it here tho
// Replaced with remove_indices_fromCmatrix above. Used in discard() methods in
// the embeddedcauset.cpp otherwise
//From https://codereview.stackexchange.com/questions/206686/removing-by-indices-several-elements-from-a-vector
template <typename INT, typename T> // INT could be int, unsigned int, char, size_t, etc...
void remove_indices(std::vector<T>& v, const std::vector<INT>& rm )
{
  // For speed, I am assuming that 'rm' is sorted
  size_t rm_index = 0;
  v.erase(
    std::remove_if(std::begin(v), std::end(v), [&](T& elem)
            {
                if (rm.size() != rm_index && &elem - &v[0] == rm[rm_index])
                {
                rm_index++;
                return true;
                }
                return false;
            }
            ), std::end(v));
}


/** @brief Gets indexes where given vector's entries have value i.
 *
 *  @param v Vector to look in. Type: std::vector<whatever>.
 *  @param x Value to look for. Type: whatever. 
 * 
 *  @return Vector of indexes. Type: vector <int>
 **/
template<typename T>
inline
std::vector <int> getIndexes(std::vector<T> v, T x)
{
    std::vector <int> sol;
    auto it = v.begin();
    while (it != v.end()) 
    {
        auto it2 = find(it+1, v.end(), x);
        if (it2 != v.end()) 
        {
        int index = it2 - v.begin();
        sol.push_back(index);
        }
        it = it2;
    }
    return sol;
}

template <typename T1, typename T2>
inline
bool contains(std::vector <T1> v, T2 x)
    {return std::find(v.begin(), v.end(), x) != v.end();
}


//_____________________________________________________________________
//
//---------------------------------------------------------------------
// VECTOR MATH FUNCTIONS
//---------------------------------------------------------------------
//_____________________________________________________________________

template <typename T1>
inline
double mymean(std::vector <T1> x)
    {return (double)myvecsum(x) / x.size();}

template <typename T1, typename F>
inline
double mymean(std::vector <T1> x, F func)
    {return (double)myvecsum(x, func) / x.size();}

template <typename T1, typename T2>
inline
double mymean(std::vector <T1> x, std::vector <T2> w = {1})
{   
    if (x.size() != w.size())
    {
        std::cout<<"Vector and weight have different sizes."<<std::endl;
        throw std::invalid_argument
        ("vector and weight have different sizes: "
        +std::to_string(x.size()) +" and " 
        +std::to_string(w.size()) + "!");
    }
    T2 sumw = myvecsum(w);
    double mean = 0;
    for (int i = 0; i < x.size(); i++)
        {mean += (double) x[i] * w[i] / sumw;}
    return mean;
}

template <typename T1, typename T2, typename F>
inline
double mymean(std::vector <T1> x, F func, std::vector <T2> w = {1})
{   
    if (w.size () == 1) return mymean(x, func);

    else if (x.size() != w.size())
    {
        std::cout<<"Vector and weight have different sizes."<<std::endl;
        throw std::invalid_argument
        ("vector and weight have different sizes: "
        +std::to_string(x.size()) +" and " 
        +std::to_string(w.size()) + "!");
    }

    T2 sumw = myvecsum(w);
    double mean = 0;
    for (int i = 0; i < x.size(); i++)
        {mean += (double) func(x[i]) * w[i];}
    return (double) mean / sumw;
}

template <typename T>
inline
double mystd(const std::vector<T>& v) {
    // Calculate the mean of the vector
    double mean = 0.0;
    for (auto x : v) {
        mean += (double)x;
    }
    mean /= v.size();

    // Calculate the sum of the squared differences from the mean
    double squaredDifferencesSum = 0.0;
    for (double x : v) {
        squaredDifferencesSum += std::pow(x - mean, 2);
    }

    // Calculate the variance and standard deviation
    double variance = squaredDifferencesSum / v.size();
    double stdDeviation = std::sqrt(variance);

    return stdDeviation;
}


inline
std::pair<double, double> combine_meass(const std::vector<int>& Ns,
                                      const std::vector<double>& mus, 
                                      const std::vector<double>& stds) 
{
    if (mus.size() != stds.size()) {
        return std::make_pair(0.0f, 0.0f); 
    }

    std::vector<double> coeffs;
    size_t M = Ns.size();
    size_t N = 0;
    for (int i = 0; i < M; i++) {N += Ns[i];}

    for (int i = 0; i < M; i++) {
        int Ni = Ns[i];
        double mui = mus[i];
        double stdi = stds[i];
        double term_i = (Ni-1)/(N-1)*stdi*stdi;
        coeffs.push_back(term_i);
        for (int j = i; j < M; j++) {
            int Nj = Ns[j];
            double muj = mus[j];
            double term_mixed = Ni*Nj/(N*(N-1)) * std::pow(mui - muj, 2);
            coeffs.push_back(term_mixed);
        }
    }
    double mu = 0.0;
    double sum_coeffs = 0.0;
    for (int i = 0; i < M; i++) {
        mu += mus[i];
    }
    mu /= M;
    for (double coeff : coeffs) {
        sum_coeffs += coeff;
    }
    double std = std::sqrt(sum_coeffs);
    return std::make_pair(mu, std);
}



/**
 * @brief Matrix multiplication A*B = C. SQUARE MATRIX PLEASE
 * 
 * @tparam T - any number type
 * @param A - matrix A
 * @param B  - matrix B
 * @return std::vector<std::vector<T>> 
 */
template <typename T>
inline
std::vector<std::vector<T>> matmul(std::vector<std::vector<T>> A,
                                   std::vector<std::vector<T>> B)
{
    // Create matrix of correct size to be filled in with results
    std::vector<std::vector<T>> C;
    int rows = A.size();
    int columns = A[0].size();
    int n = rows;
    if (rows != columns) {
        print("Want square matrix please.");
        throw std::runtime_error("");
    }
    if ((B.size()!=n) || (B[0].size()!=n)) {
        print("A and B matrix should be of equal sizes!");
        throw std::runtime_error("");

    }
    C.resize(rows, std::vector<T>(columns));

    //#pragma omp parallel for schedule(dynamic)
    //#pragma omp parallel for private(i,j,k) shared(A,B,C)    
    for(int i = 0; i < n; i++) {        
        for(int k = 0; k < n; k++) {       
            for(int j = 0; j < n; j++) {                
                C[i][j] += A[i][k] * B[k][j];  
            }
        }
    }
    return C;
}

/**
 * @brief Sum of all elements in the matrix (assuming a rectangular matrix!)
 * 
 * @tparam T - any real number type
 * @param M matrix to be summed over
 * @return double 
 */
template <typename T>
inline
double sumMatrix(std::vector<std::vector<T>> M)
{
    double sum = 0;
    int N = M.size();
    int K = M[0].size();
    for (int i=0; i<N; i++){
        for (int j=0; j<K; j++) {
            sum = sum + M[i][j];
        }
    }
    return sum;
}

// CTRL+K+C to comment out multiple lines of code
// CTRL+K+U to uncomment out multiple lines of code

// //_____________________________________________________________________
// //
// //---------------------------------------------------------------------
// // ADVANCED FUNCTIONS
// //---------------------------------------------------------------------
// //_____________________________________________________________________
// template <typename T1, typename T2>
// inline
// vector <T1> getAwhereB(vector<T1> v1, vector<T2> v2, T2 x)
// /** @brief Gets values of vector v1 at indexes where vector v2 = x.
//  *
//  *  @param v1 Vector to take values from. Type: std::vector<whatever>.
//  *  @param v2 Vector where to look for. Type: std::vector<whatever2>. 
//  *  @param x  Value to look for. Type: whatever2. 
//  * 
//  *  @return Values of vector v1 at indexes where vector v2 = x.
//  **/
// {   
//     if (v1.size() != v2.size())
//     {
//         throw std::invalid_argument( "Vectors Have Different Sizes" );
//     } 
//     vector <T1> sol;
//     auto begin = v2.begin();
//     auto end = v2.end();
//     auto it = begin-1;
//     while (it != end) 
//     {
//         auto it2 = find(it+1, end, x);
//         if (it2 != end) 
//         {
//         int index = it2 - begin;
//         sol.push_back(v1[index]);
//         }
//         it = it2;
//     }
//     return sol;
// }

// template <typename T1, typename T2>
// inline
// vector<vector<T1>> getAwhereB(vector<vector<T1>> V, vector<T2> u, T2 x)
// /** @brief Gets values of vectors v0=V[0], v1=V[1], etc... at indexes 
//  *         where vector U = x.
//  *
//  *  @param V  Vector of many vectors to take values from. 
//  *            Type: std::vector< vector<whatever> >.
//  *  @param U Vector where to look for. Type: std::vector<whatever2>. 
//  *  @param x  Value to look for. Type: whatever2. 
//  * 
//  *  @return Vector of vectors with values corrsponding to indexes where
//  *          vector U = x.
//  **/
// {   
//     for (int i = 0; i < V.size(); i++)
//     {
//         if (V[i].size() != u.size())
//         {
//             throw std::invalid_argument
//             ( "Vector " + to_string(i) + " in V and Vector U have Different Sizes" );
//         } 
//     }
//     vector< vector<T1> > sol (V.size()); //one vector per each vector in V
//     auto begin = u.begin();
//     auto end = u.end();
//     auto it = begin-1;
//     while (it != end) 
//     {
//         auto it2 = find(it+1, end, x);
//         if (it2 != end) 
//         {
//         int index = it2 - begin;
//         for (int i = 0; i < V.size(); i++)
//             {
//                 sol[i].push_back(V[i][index]);
//             }
//         }
//         it = it2;
//     }
//     return sol;
// }

// template <typename T1, typename T2, typename F>
// vector <T1> getAwhereB(vector<T1> v1, vector<T2> v2, F f)
// /** @brief Gets values of vector v1 at indexes where v2 satisfies f.
//  *
//  *  @param v1 Vector to take values from. Type: std::vector<whatever>.
//  *  @param v2 Vector where to look for. Type: std::vector<whatever2>. 
//  *  @param f  Function taking as argument whatever2 type and 
//  *            returning bool. 
//  * 
//  *  @return Values of vector v1 at indexes where vector v2 satisfies f.
//  **/
// {   
//     if (v1.size() != v2.size())
//     {
//         throw std::invalid_argument( "Vectors Have Different Sizes" );
//     } 
//     vector <T1> sol;
//     auto begin = v2.begin();
//     auto end = v2.end();
//     auto it = begin-1;
//     for (int ind = 0; ind < v2.size(); ind++) 
//         {if ( f(v2[ind]) == 1 ) sol.push_back(v1[ind]);}
//     return sol;
// }

// template <typename T1, typename T2, typename F>
// vector <vector<T1>> getAwhereB(vector<vector<T1>> v1, vector<T2> v2, F f)
// /** @brief Gets values of vector v1 at indexes where v2 satisfies f.
//  *
//  *  @param v1 Vector to take values from. Type: std::vector<whatever>.
//  *  @param v2 Vector where to look for. Type: std::vector<whatever2>. 
//  *  @param f  Function taking as argument whatever2 type and 
//  *            returning bool. 
//  * 
//  *  @return Values of vector v1 at indexes where vector v2 satisfies f.
//  **/
// {   
//     if (v1.size() != v2.size())
//     {
//         throw std::invalid_argument( "Vectors Have Different Sizes" );
//     } 
//     vector <T1> sol;
//     auto begin = v2.begin();
//     auto end = v2.end();
//     auto it = begin-1;
//     for (int ind = 0; ind < v2.size(); ind++) 
//         {if ( f(v2[ind]) == 1 ) sol.push_back(v1[ind]);}
//     return sol;
// }


// template <typename T1, typename T2>
// vector <T1> sortAwithB(vector <T1> A, vector <T2> B, bool reverse = false)
// {
//     struct mypair
//     {
//         T1 a;
//         T2 b;
//     };

//     vector <mypair> zipped (A.size());
//     for (int i = 0; i<A.size(); i++)
//     {
//         zipped[i] = {A[i], B[i]};
//     }

//     if (not reverse)
//     sort(zipped.begin(), zipped.end(), 
//          [](const auto& i, const auto& j){ return i.b < j.b; } );
//     else
//     sort(zipped.begin(), zipped.end(), 
//          [](const auto& i, const auto& j){ return i.b > j.b; } );
    
//     vector <T1> sortA (A.size());
//     for (int i = 0; i<A.size(); i++)
//     {
//         sortA[i] = zipped[i].a;
//     }
    
//     return sortA;
// }


// template <typename T1>
// vector <T1> sum_portions(vector <T1> v, int n)
// /** @brief Sums first, second, ..., nth portion of vector 
//  *         with each other.
//  *
//  *  @param v Vector. Type vector <whatever>.
//  *  @param n Number of portions into which divide v, then sum
//  *           on each other. Type int.

//  *  @return Vector of v.size()/n elements. 
//  **/
// {
//     if (v.size()%n) 
//     { 
//         throw std::invalid_argument("Vector's size not multiple of n");
//     }
//     int m = v.size()/n;
//     vector <T1> sol (m); 
 
//     for (int i = 0; i < v.size(); i++) 
//     {
//         int index = i%m;
//         sol[index] += v[i];
//     }
//     return sol;
// }


// template <typename T1, typename T2>
// inline
// double minAbsDiff (vector<T1> v1, vector<T2> v2)
// /** @brief Minimum absolute difference between elements of v1 and v2.
//  *
//  *  @param v1 vector 1 of numbers
//  *  @param v2 vector 2 of numbers
//  * 
//  *  @return Minimum absolute difference between elements of v1 and v2.
//  *          Type: double.
//  **/
// {
//     double min = INFINITY;  
//     for (int i=0; i<v1.size(); i++)
//     {
//         for (int j=0; j<v2.size(); j++)
//         {
//             double diff = abs( v1[i] - v2[j] );
//             min = (diff < min) ? diff : min;
//         }
//     }
//     return min;
// }


// //https://gist.github.com/lorenzoriano/5414671
// template <typename T>
// std::vector<T> linspace(T a, T b, size_t N) {
//     T h = (b - a) / static_cast<T>(N-1);
//     std::vector<T> xs(N);
//     typename std::vector<T>::iterator x;
//     T val;
//     for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
//         *x = val;
//     return xs;
// }




/*
template <typename T>
inline
T cumulative (vector <T> v)
{
    vector <double> cum (N);
    T tot = myvecsum(v);
    v[0] /= tot;
    partial_sum(v.begin(), v.end(), cum.begin(), 
               [](double a, double b)
                {return (double)a + b;});
    return cum;
}

template <typename T>
inline
double cumulative (vector <T> v, bool norm = true)
{
    vector <double> cum (N);
    T tot = myvecsum(v);
    v[0] /= tot;
    if (norm) auto lambda = [tot](double a, double b)
                            {b/=tot; return (double)a + b;};
    else      auto lambda = [](double a, double b)
                            {return (double)a + b;};
    partial_sum(v.begin(), v.end(), cum.begin(), lambda);
    return cum;
}*/

#endif /* MYCFUNCTIONS_H */