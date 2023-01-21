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

#include "../scripts_cpp/causets_cpp/causet.h"

// Use the "bisection" root finder method
//using boost::math::tools::bisect;
using std::cout;
using std::endl;
using std::vector;

template <typename F>
double bisection (F f, double xleft, double xright, double epsilon = 1e-8, 
                  int Nstop = 1e3)
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


double MM_drelation(double d)
{
    double a = std::tgamma(d+1);
    double b = std::tgamma(d/2);
    double c = 4* std::tgamma(3*d/2);
    return a*b/c;
}

double func(double x)
{
  return x*x*x-8;
}


int main()
{
  double from = 0.9;  
  double to = 7.1;


  // double soln = 2;
  // double result = bisection(func, from, to); 
  // cout << "Result " << result << "  Correct result: "<< soln << endl;
  // cout << "The difference is " << result - soln << endl;
  // cout << endl;

  //vector<double> dims = {0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5};
  //for (auto x: dims)
  //{
  //  double result =  MM_drelation(x);
  //  cout << "Dim =" << x << " Func-cpp = " << result << endl;
  //}



vector<double> frs = {0.5, 0.426169, 0.359442, 0.300726,
      0.25, 0.206757, 0.170263, 0.139708,
      0.11428571428, 0.093243, 0.0759003, 0.0616587,
      0.05, 0.0404813, 0.032728, 0.0264256,0.02131202131};


vector<double> sols = {1, 1.25, 1.5, 1.75, 2,
                       2.25, 2.5, 2.75, 3, 3.25,
                       3.5, 3.75, 4, 4.25,4.5,4.75,5};

}

    

