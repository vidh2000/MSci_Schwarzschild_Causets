#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <set>
#include <stack>
#include <stdio.h>
#include <stdexcept>
#include <string>
#include <vector>
#include <chrono>
#include <unordered_set>
#include <random>
#include <
// Use the "bisection" root finder method
using boost::math::tools::bisect;


struct TerminationCondition  {
  bool operator() (double min, double max)
  {
    return abs(min - max) <= 0.000001;
  }
};

double FunctionToApproximate(double x)
{
    return (x*x*x - 8); //root at x=2
}


int main(){
// The solution must lie in the interval [from, to], additionally f(from) <= 0 && f(to) >= 0
double from = 0;  
double to = 6;
std::pair<double, double> result = bisect(FunctionToApproximate(),
                                          from, to, TerminationCondition());
double root = (result.first + result.second) / 2;  

std::cout << "Result: x =" << root << ". Correct result: x = 2" << std::endl;

}

