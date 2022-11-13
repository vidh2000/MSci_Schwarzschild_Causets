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

#include "boost/math/tools/roots.hpp"

double poly (double x){
    return x*x*x - 8;
};

struct TermCond  {
  bool operator() (double min, double max)  {
    return abs(min - max) <= 0.000001;
  }
};


int main()
{
    auto r = boost::math::tools::bisect(poly, 0, 8, TermCond());
    std::cout<<r.first<<", "<<r.second<<std::endl;
    return 0;
}