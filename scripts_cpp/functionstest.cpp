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

#include "functions.h"
#include "vecfunctions.h"
using namespace std;

int main()
{
    // std::set<int> a = {1,2,3,4};
    // std::set<int> b = {3,4,5};
    // a = set_add(a,b);
    // print_set(a);

    vector<double> v = linspace(0., 100., 100);
    print_vector(v);    
}