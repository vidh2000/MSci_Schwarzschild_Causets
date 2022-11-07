#include <algorithm>
#include <boost/math/tools/roots.hpp>
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
#include <boost\math\tools\roots.hpp> 


// Use the "bisection" root finder method
using boost::math::tools::bisect;
using std::vector;
using std::cout;
using std::endl;


struct TerminationCondition  
{
  bool operator() (double min, double max)
  {
    return abs(min - max) <= 0.000001;
  }
};

double FunctionToApproximate(double x)
{
  return (x*x*x - 8); //root at x=2
}

double MM_drelation(double d)
{
    double a = std::tgamma(d+1);
    double b = std::tgamma(d/2);
    double c = 4* std::tgamma(3*d/2);
    return a*b/c;
}

double tol(double min, double max)
  {return max - min > 1e-8;}


int main()
{
// The solution must lie in the interval [from, to], additionally f(from) <= 0 && f(to) >= 0
double from = 0.75;  
double to = 6;

vector<double> frs = {1/2, 0.426169, 0.359442, 0.300726, 1/4, 0.206757, 0.170263, 0.139708, 4/35, 0.093243, 0.0759003, 0.0616587, 1/20, 0.0404813, 0.032728, 0.0264256, 64/3003};
vector<double> sols = {1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.25,4.5,4.75,5};
double fr_i;
double sol_i;
for (int i = 0; i<frs.size(); i++)
{
  fr_i = frs[i];
  auto dim_solve = [fr_i](double x){return MM_drelation(x) - fr_i;};
  std::pair<double, double> bisect_pair = bisect (dim_solve, from, to, tol);
  double result = (bisect_pair.second - bisect_pair.first)/2;
  cout << "Result " << result << "  Correct result: "<< sols[i] << endl;
  cout << " The difference is " << result - sols[i] << endl;
  //cout << endl;

}

} /*end main*/




    

