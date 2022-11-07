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
        std::string separator = " , ";
        std::cout << e << separator;
    }
    std::string endstr = " }";
    std::cout << endstr;
    std::cout << std::endl;
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


//for set containing int only here
template <typename obj> 
inline
std::set<obj> set_diff(std::set<obj> s1, std::set<obj> s2)
    /*
    Return set difference i.e events which are only in s1;
    RETURN = s1-s2 (where s2 can have other elements as well)
    */
{
    std::set<obj> result;
    std::set_difference(s1.begin(), s1.end(), s2.begin(), s2.end(),
        std::inserter(result, result.end()));
    return result;
}

template <typename obj> 
inline
std::set<obj> set_union(std::set<obj> s1, std::set<obj> s2)
    /**
    * @brief Returns union of two sets s1 and s2.
    */
{
    s1.insert(s2.begin(), s2.end());
    return s1;
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


#endif /* FUNCTIONS_H */