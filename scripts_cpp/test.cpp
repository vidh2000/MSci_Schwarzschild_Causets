#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <set>
#include <stack>
#include <stdio.h>
#include <stdexcept>
#include <string>
#include <vector>
#include <chrono>
#include <unordered_set>
#include <chrono>
#include <typeinfo>


#include "causets_cpp/sprinkledcauset.h"
#include "causets_cpp/shapes.h"
#include "causets_cpp/spacetimes.h"

#include "causets_cpp/functions.h"
#include "causets_cpp/vecfunctions.h"



int main(){

    std::unordered_set<int> a = {1,2,5,13,8};
    std::unordered_set<int> b = {2,4,5,6,8,10,13};
    print_set(a);
    print_set(b);
    std::unordered_set<int> intersection = set_intersection(a,b);
    print_set(intersection);

    return 0;
}