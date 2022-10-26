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

inline
void print_set(std::set<int> set)
{
    for (int const& e : set)
    {
        std::cout << e << ' ';
    }
    std::cout << std::endl;
}

inline
bool set_contains(int element, std::set<int> s)
{
    bool is_in = s.find(element) != s.end();
    return is_in;
}


//for set containing int only here
inline
std::set<int> set_diff(std::set<int> s1, std::set<int> s2)
    /*
    Return set difference i.e events which are only in s1;
    RETURN = s1-s2 (where s2 can have other elements as well)
    */
{
    std::set<int> result;
    std::set_difference(s1.begin(), s1.end(), s2.begin(), s2.end(),
        std::inserter(result, result.end()));
    return result;
}

inline
std::set<int> set_add(std::set<int> s1, std::set<int> s2)
    /**
    * @brief Returns union of two sets s1 and s2.
    */
{
    s1.insert(s2.begin(), s2.end());
    return s1;
}


#endif /* FUNCTIONS_H */