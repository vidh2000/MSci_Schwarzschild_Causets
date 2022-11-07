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
#include <stdlib.h>     /* abs */
#include <stdexcept>
#include <string>
#include <vector>


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
    std::uniform_real_distribution<> dis(0,1.0);

    for(int i=0; i<rows; i++)
    {
        matrix[i].resize(cols);
        for (int j=0; j<cols; j++){
            matrix[i][j] = dis(gen);
        }
    }
    return matrix;
}


/**
 * @brief Generate "size" DIFFERENT random ints between 0 and N-1.
 */
inline
std::vector<int> distinct_randint1(int size, int N, int seed)
{
    std::vector<int> result(size);
    if (!seed)
    {
        std::random_device rd;
        seed = rd();
    }  
    std::mt19937 gen(seed);
    std::uniform_real_distribution<> dis(0,1.0);
    
    for(int i = 0; i < size; ++i)
    {
        int r;
        while(std::find(result.begin(), result.end(), r) != result.end())
        {
            r = dis(gen) * N;
            result[i] = r;   
        }
    return result;
    }
}

/**
 * @brief Generate "size" DIFFERENT random ints between 0 and N-1.
 */
inline
std::vector<int> distinct_randint2(int size, int N, int seed)
{
    std::vector<int> result(N); 
    if (!seed)
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
std::vector<int> distinct_randint(int size, int N, int seed)
{
    if (size < N/2){
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
        if (brackets) {std::cout<<"\n{ ";}
        for(int index = 0; index < size ; index++)
        {   
            if (index == 0) {std::cout<<vec[index];}
            else {std::cout<<sep<<vec[index];}
        }
        if (brackets) {std::cout<<" }";}
        std::cout<<std::endl;  
}


template <typename T>
inline
void print_vector(std::vector<std::vector<T>> vec, 
                  std::string sep=" , ",
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
        if (brackets) {std::cout<<"\n{ ";}
        for(int index = 0; index < size ; index++)
        {
            std::vector<T> vind = vec[index];
            std::cout<<"{";
            for (int j = 0; j < vind.size(); j++)
            {
                if (j == 0) {std::cout<<vind[0];}
                else {std::cout<<sep<<vind[j];}
            }
            std::cout<<"}";

            if (index != size -1) {std::cout<<sep;}
        }
        if (brackets) {std::cout<<" }";}
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
        throw std::invalid_argument
        ("vector and weight have different sizes: "
        +to_string(x.size()) +" and " 
        +to_string(w.size()) + "!");
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
        throw std::invalid_argument
        ("vector and weight have different sizes: "
        +to_string(x.size()) +" and " 
        +to_string(w.size()) + "!");
    }

    T2 sumw = myvecsum(w);
    double mean = 0;
    for (int i = 0; i < x.size(); i++)
        {mean += (double) func(x[i]) * w[i];}
    return (double) mean / sumw;
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