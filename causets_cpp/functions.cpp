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

void print_set(std::set<int> set)
{
    std::cout << "{ ";
    for (int const& e : set)
    {
        std::cout << e << ' ';
    }
    std::cout << std::endl;
}
bool set_contains(int element, std::set<int> s)
{
    bool is_in = s.find(element) != s.end();
    return is_in;
}


//for set containing int only here
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

std::set<int> set_union(std::set<int> s1, std::set<int> s2)
    /*
    Returns union of two sets
    RETURN = s1 U s2
    */
    {
    s1.insert(s2.begin(), s2.end());
    return s1;
    }

/*
int main(){
std::set<int> a = {1,2,3,4};
std::set<int> b = {3,4,5};
a = set_add(a,b);
print_set(a);

}
*/