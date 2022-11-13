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
#include "boost/numeric/odeint.hpp"

using namespace boost::numeric::odeint;
using std::cout;
using std::endl;

double poly (double x){
    return x*x*x - 8;
};

struct TermCond  {
  bool operator() (double min, double max)  {
    return abs(min - max) <= 0.000001;
  }
};

double BH_dvarphi_du_plus (double u, double M, double c2)
{
  return std::pow(2*M*u*u*u - u*u + c2, -0.5);
};
double BH_dvarphi_du_minus (double u, double M, double c2)
{
  return - std::pow(2*M*u*u*u - u*u + c2, -0.5);
};
double BH_dt_du_plus (double u, double M, double c)
{
  return (c* std::pow(2*M*u*u*u - u*u + c*c, -0.5) - 2*M*u)/
                       (u*u - 2*M*u*u*u);
};
double BH_dt_du_minus (double u, double M, double c)
{
  return (-c* std::pow(2*M*u*u*u - u*u + c*c, -0.5) - 2*M*u)/
                       (u*u - 2*M*u*u*u);
};
double c2_solver (double u1, double u2, double varphi2, double M)
{
  double D1 = u1*u1 - 2*M*u1*u1*u1;
  double D2 = u2*u2 - 2*M*u2*u2*u2;
  double Dmax = (D1>D2)? D1 : D2;
  double c2min = (Dmax>0)? Dmax : 0;
  if (u2>u1 <0)
  {
    // USE PLUS
    double c2max = (u2+(u2-u1)/varphi2)*(u2+(u2-u1)/varphi2); 
    // SOLVE WITH BISECTION
  }
  else
  {
    //USE MINUS
    double c2max = (u1+(u1-u2)/varphi2)*(u1+(u1-u2)/varphi2); 
    // SOLVE WITH BISECTION
  }
};


int main()
{

    cout<<"\n============ TESTING BISECTION =============\n";
    auto r = boost::math::tools::bisect(poly, 0, 8, TermCond());
    cout<<"Boost's solution to x^3-8 is "<<(r.first+r.second)/2<<endl;

    cout<<"\n============ TESTING DIFF EQS =============\n";
    double u1 = 1./2.;
    double u2 = 1./3.;
    double deltau = u1-u2;
    double varphi1 = 0;
    double varphi2 = 1.;
    size_t steps = integrate(BH_dvarphi_du_plus, varphi1, u1, u2, deltau/50);
    return 0;
}