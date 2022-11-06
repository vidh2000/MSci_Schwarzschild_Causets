#include <algorithm>
//#include <boost/>
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

#include "../scripts_cpp/causets_cpp/causet.h"

// Use the "bisection" root finder method
//using boost::math::tools::bisect;


template <typename F>
double bisection (F f, double xleft, double xright, double epsilon = 1e-5, 
                  int Nstop = 1e3)
{
  double fl = f(xleft);
  double fr = f(xright);
  double xnew = 0;

  int counter = 0;
  double e = 1;
  while (e > epsilon && counter < Nstop)
  {
      xnew = (xleft + xright) /2;
      double fnew = f(xnew);

      if (fnew * fl >= 0)
      {
        e = xleft-xnew;
        e = (e>0)? e : -e;
        xleft = xnew*1;
      }
      else
      {
        e = xright-xnew;
        e = (e>0)? e : -e;
        xright = xnew*1;
      }
      counter ++;
  }
  return xnew;
}

double FunctionToApproximate(double x)
{
  return (x*x*x - 8); //root at x=2
}


int main()
{
// The solution must lie in the interval [from, to], additionally f(from) <= 0 && f(to) >= 0
std::cout<<abs(0.3)<<"\n";
double from = 0;  
double to = 6;
double result = bisection(FunctionToApproximate,
                                              from, to); 

std::cout << "Result: x =" << result << ". Correct result: x = 2" << std::endl;
}

    

