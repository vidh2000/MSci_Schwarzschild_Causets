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


#endif /* FUNCTIONS_H */