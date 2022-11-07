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
using std::cout;
using std::endl;

template <typename F>
double bisection (F f, double xleft, double xright, double epsilon = 1e-4, 
                  int Nstop = 1e4)
{
  double fl = f(xleft);
  double fr = f(xright);
  double xnew = 0;
  double fnew = 0;

  int counter = 0;
  double e = 1;
  while (e > epsilon && counter < Nstop)
  {
      xnew = (xleft + xright) /2;
      fnew = f(xnew);

      if (fnew * fl >= 0)
      {
        e = xnew-xleft;
        xleft = xnew*1;
      }
      else
      {
        e = xright-xnew;
        xright = xnew*1;
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



int main()
{
double from = 0.5;  
double to = 5;
vector<double> frs = {3.1415/8, 0.387412, 0.461414, 0.638338, 1, 1.74226, 
                      3.33083, 6.91527, (315 *3.1415)/64, 36.9854, 94.1055, 
                      253.496, 720, 2148.66};
vector<double> sols = {1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 
                       3, 3.25, 3.5, 3.75, 4, 4.25};

double fr_i;
double sol_i;
for (int i = 0; i<frs.size(); i++)
{
  fr_i = frs[i];
  auto dim_solve = [fr_i](double x){return MM_drelation(x) - fr_i;};
  double result = bisection(dim_solve, from, to); 
  cout << "Result " << result << "  Correct result: "<< sols[i] << endl;
  cout << "The difference is " << result - sols[i] << endl;
  cout << endl;
}
}

    

