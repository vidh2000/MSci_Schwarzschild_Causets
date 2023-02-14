#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <fstream>
#include <stack>
#include <string>
#include <stdio.h>
#include <vector>
#include <set>
#include <unordered_set>

/**
 * @brief           Root finder in some interval [x_left,x_right]
 *                  using the "bisection" iterative method.
 *                  Requires only 1 root to be in the interval for
 *                  the function to work correctly and the function
 *                  NEEDS to be monotone in the provided interval.
 * 
 * @tparam F :      function/method returning a double
 * @param f:        the function name whose roots we're searching for
 * @param xleft:    the left boundary of search
 * @param xright:   the right boundary of search
 * @param epsilon:  the allowed uncertainty (error) in the result 
 * @param Nstop:    the upper limit of iterations
 * @return double 
 */
template <typename F>
inline
double bisection (F f, double xleft, double xright, double epsilon = 1e-8, 
                  int Nstop = 1000)
{
  double fl = f(xleft);
  double fr = f(xright);
  double xnew;
  double fnew;

  int counter = 0;
  double e = 1;
  while (e > epsilon && counter < Nstop)
  {
      xnew = (xleft + xright) /2;
      fnew = f(xnew);
      if (fnew * fl > 0)
      {
        e = xnew-xleft;
        xleft = xnew*1.0;
      }
      else
      {
        e = xright-xnew;
        xright = xnew*1.0;
      }
      counter ++;
  }
  return xnew;
}

inline
double MM_drelation(double d)
{
    double a = std::tgamma(d+1);
    double b = std::tgamma(d/2);
    double c = 4* std::tgamma(3*d/2);
    return a*b/c;
}

template <typename obj>  
inline
void print_set(std::set<obj> set)
{
    
    std::string beginstr = "{ ";
    std::cout << beginstr;
    for (obj e : set)
    {
        std::string separator = " ";
        std::cout << e << separator;
    }
    std::string endstr = "}";
    std::cout << endstr;
    std::cout << std::endl;
}

template <typename obj>  
inline
void print(std::set<obj> set)
{
    
    std::string beginstr = "{ ";
    std::cout << beginstr;
    for (obj e : set)
    {
        std::string separator = " ";
        std::cout << e << separator;
    }
    std::string endstr = "}";
    std::cout << endstr;
    std::cout << std::endl;
}


template <typename obj>  
inline
void print_set(std::unordered_set<obj> set)
{
    
    std::string beginstr = "{ ";
    std::cout << beginstr;
    for (obj e : set)
    {
        std::string separator = " ";
        std::cout << e << separator;
    }
    std::string endstr = "}";
    std::cout << endstr;
    std::cout << std::endl;
}

template <typename obj>  
inline
void print(std::unordered_set<obj> set)
{
    
    std::string beginstr = "{ ";
    std::cout << beginstr;
    for (obj e : set)
    {
        std::string separator = " ";
        std::cout << e << separator;
    }
    std::string endstr = "}";
    std::cout << endstr;
    std::cout << std::endl;
}


inline
void print(std::string str)
{
    std::cout << str << std::endl;
}

template <typename obj>  
inline
void print(std::vector<std::unordered_set<obj>> vec)
{
    bool brackets = true;
    int size = vec.size();
    if (brackets) {std::cout<<"{";}
    for(int index = 0; index < size ; index++)
    {
        std::unordered_set<obj> set = vec[index];
        if (index != 0){
            std::cout<<" {";
        }
        else {
            std::cout << "{";
        }
        for (auto e : set)
        {
            std::string separator = " ";
            std::cout << e << separator;
        }
        std::cout << "}";

        if (index != size -1) {std::cout<< ", " << std::endl;}
    }
    if (brackets) {
        std::cout<<"}";
    }
    std::cout<<std::endl;
}


template <typename obj>  
inline
bool set_contains(obj element, std::set<obj> s)
{
    bool is_in = s.find(element) != s.end();
    return is_in;
}

template <typename obj>
inline
bool set_contains(obj element, std::unordered_set<obj> s)
{
    bool is_in = s.find(element) != s.end();
    return is_in;
}



template <typename E> 
inline
std::set<E> set_diff(std::set<E> s1, std::set<E> s2)
    /*
    Return set difference i.e events which are only in s1;
    RETURN = s1-s2 (where s2 can have other elements as well)
    */
{
    std::set<E> result;
    std::set_difference(s1.begin(), s1.end(), s2.begin(), s2.end(),
        std::inserter(result, result.end()));
    return result;
}


template <typename E> 
inline
std::unordered_set<E> set_diff(std::unordered_set<E> s1,
                        std::unordered_set<E> s2)
    /*
    Return set difference i.e events which are only in s1;
    RETURN = s1-s2 (where s2 can have other elements as well)
    */
{
    for (const auto& elem : s2) {
    s1.erase(elem);
    }
    return s1;
}



template <typename SET> 
inline
SET set_union(SET s1, SET s2)
    /**
    * @brief Returns union of two sets s1 and s2.
    */
{
    s1.insert(s2.begin(), s2.end());
    return s1;
}

/**
 * @brief Finds intersection of two sets s1, s2.
 *        -- All the elements that belong to both
 *           sets simultaneously.
 * 
 * @param s1: set 1
 * @param s2: set 2
 * @return Intersection of s1, s2 of type SET
 */
template <typename obj>
inline
std::set<obj> set_intersection(std::set<obj> s1, std::set<obj> s2)
{
    std::set<obj> result;
    std::set_intersection(s1.begin(), s1.end(),s2.begin(), s2.end(),
                          std::inserter(result,result.end()));
    return result;
}

/**
 * @brief Finds intersection of two sets s1, s2.
 *        -- All the elements that belong to both
 *           sets simultaneously.
 * 
 * @param s1: set 1
 * @param s2: set 2
 * @return Intersection of s1, s2 of type SET
 */
template <typename obj>
inline
std::unordered_set<obj> set_intersection(std::unordered_set<obj> s1,
                                         std::unordered_set<obj> s2)
{
    // Require sorted object...

    std::unordered_set<obj> result;
    for (int element : s1)
    {
        if (s2.count(element) > 0) {
        result.insert(element);
        }
    }
    return result;
}

// Function to find the maximum element
template <typename num> 
inline
num setmax(std::set<num> my_set)
{
 
    // Get the maximum element
    int max_element;
    if (!my_set.empty())
        max_element = *(my_set.rbegin());
 
    // return the maximum element
    return max_element;
}
 
// Function to find the minimum element
template <typename num> 
inline
num setmin(std::set<num> my_set)
{
 
    // Get the minimum element
    int min_element;
    if (!my_set.empty())
        min_element = *my_set.begin();
 
    // return the minimum element
    return min_element;
}

inline
void discard_from_set(std::unordered_set<int> &myset, int label)
{
    int N = myset.size();
    std::unordered_set<int> buffer;
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
inline
void discard_from_set(std::unordered_set<m> &myset, std::vector<int> labels)
{
    std::sort(labels.begin(),labels.end());
    int N = myset.size();
    std::unordered_set<m> buffer;
    int startpoint = 0;
    for (m j : myset) //not ordered
    {
        if (j>labels[-1])
            {buffer.insert(j-labels.size());}
        else
        {
            for (int s = 0; s<labels.size(); s++)
                {
                    if (labels[s] == j)
                        {break;}
                    else if (labels[s] > j)
                        {buffer.insert(j-s); break;}
                }
        }
    }
    myset = buffer;
}


/**
 * @brief Replaces the values in the vector of (unordered) sets
 *        with indices corresponding to the position where those
 *        values appear in the "interval" vector
 * @param sets : i.e futures/pasts. Vector of sets which contain ints etc.
 * @param interval: vector which contains the target values and which will 
 *                  be used to get the indices corresponding to their loc in it
 */
template<typename T>
inline
void replace_indices(std::vector<std::unordered_set<T>> &sets,
                std::vector<T> interval)
{
    print("Inside 'replace_indices' method in functions.h");
    // Make sure values in the interval are ordered
    std::sort(interval.begin(),interval.end());
    std::vector<std::unordered_set<T>> new_sets;
    print("Updated sets:");
    for (int i = 0; i<sets.size(); i++)
    {
        if (std::find(interval.begin(),interval.end(),i) == interval.end()) {
            continue;
        }
    
        std::unordered_set<T> new_set = {};
        for (const auto& value : sets[i])
        {
            auto it = std::find(interval.begin(), interval.end(), value);
            if (it != interval.end()) {
                new_set.insert(std::distance(interval.begin(), it));
            }
        }
        new_sets.push_back(new_set);
        print(new_set);
          
    }
    sets = new_sets;
} 


/**
 * @brief Replaces the values in the set
 *        with indices corresponding to the position where those
 *        values appear in the "interval" vector
 * @param sets : i.e future/past. Set which contain ints etc.
 * @param interval: vector which contains the target values and which will 
 *                  be used to get the indices corresponding to their loc in it
 */
template<typename T>
inline
void replace_indices(std::unordered_set<T> &set,
                std::vector<T> interval)
{

    // Make sure values in the interval are ordered
    std::sort(interval.begin(),interval.end());
    
    std::unordered_set<T> new_set = {};
    for (const auto& value : set)
    {
        auto it = std::find(interval.begin(), interval.end(), value);
        if (it != interval.end()) {
            new_set.insert(std::distance(interval.begin(), it));
        }
    }
    set = new_set;
}
         




/**
 * @brief The same as np.arange
 * 
 * @tparam T - any value integer, float...
 * @param start Smallest value
 * @param stop  Up to where the interval should be
 * @param step  Step between values in the interval [start,stop)
 * @return std::vector<T> 
 */
template<typename T>
std::vector<T> arange(T start, T stop, T step = 1) {
    std::vector<T> values;
    for (T value = start; value < stop; value += step)
        values.push_back(value);
    return values;
}




#endif /* FUNCTIONS_H */