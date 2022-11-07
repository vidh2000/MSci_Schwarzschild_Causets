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

using std::vector;
using std::set;
using std::unordered_set;


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


void discard_from_set(unordered_set<int> &myset, int label)
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
void discard_from_set(unordered_set<m> &myset, vector<int> labels)
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


#endif /* FUNCTIONS_H */