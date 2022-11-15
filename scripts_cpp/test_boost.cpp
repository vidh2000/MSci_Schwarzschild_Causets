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

#include "functions.h"
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

void BH_dvarphi_du (double& dpdu, double u, double M, double c2)
{
  dpdu = std::pow(2*M*u*u*u - u*u + c2, -0.5);
};
double BH_int_dvarphi_du(double u1, double u2, double M, double c2)
{
    auto BH_dvarphi_du_forint = [M, c2]
                            (const double& varphi, double& dpdu, const double u)
                            {dpdu = std::pow(2*M*u*u*u - u*u + c2, -0.5);};
    double varphi = 0;
    if (u2>=u1)
    {
        boost::numeric::odeint::integrate(BH_dvarphi_du_forint, varphi, 
                                          u1, u2, (u2-u1)/20);
    }
    else
    {
        boost::numeric::odeint::integrate(BH_dvarphi_du_forint, varphi, 
                                          u2, u1, (u1-u2)/20);
    }
    return varphi;
};
void BH_dt_du_plus (double&dtdu, double u, double M, double c)
{
  dtdu = (c* std::pow(2*M*u*u*u - u*u + c*c, -0.5) - 2*M*u)/
                        (u*u - 2*M*u*u*u);
};
void BH_dt_du_minus (double&dtdu, double u, double M, double c)
{
  dtdu = (-c* std::pow(2*M*u*u*u - u*u + c*c, -0.5) - 2*M*u)/
                        (u*u - 2*M*u*u*u);
};
double BH_int_dt_du (double u1, double u2, double M, double c)
{
    double t = 0;
    if (u2>=u1)
    {
        auto BH_dt_du_forint_plus = [M, c]
                              (const double& t, double& dtdu, const double u)
                              {dtdu = (c* std::pow(2*M*u*u*u - u*u + c*c, -0.5) 
                                        - 2*M*u)/
                                      (u*u - 2*M*u*u*u);
                              };
        boost::numeric::odeint::integrate(BH_dt_du_forint_plus, t, 
                                          u1, u2, (u2-u1)/20);
    }
    else
    {
        auto BH_dt_du_forint_minus = [M, c]
                              (const double& t, double& dtdu, const double u)
                              {dtdu = (-c*std::pow(2*M*u*u*u - u*u + c*c, -0.5) 
                                        - 2*M*u)/
                                      (u*u - 2*M*u*u*u);
                              };
        boost::numeric::odeint::integrate(BH_dt_du_forint_minus, t, 
                                          u1, u2, (u1-u2)/20);
    }
    return t;
};

double BH_c2_solver (double u1, double u2, double varphi2, double M)
{
  //SET MINIMUM BOUND
  double D1 = u1*u1 - 2*M*u1*u1*u1;
  double D2 = u2*u2 - 2*M*u2*u2*u2;
  double Dmax = (D1>D2)? D1 : D2;
  double c2min = (Dmax>0)? Dmax +1e-4 : 0;
  // U2>=U1 -> USE PLUS
  if (u2>=u1) 
  {
    double deltaphi_max = BH_int_dvarphi_du(u1, u2, M, c2min);
    if (deltaphi_max <= 0)
      {return std::sqrt(c2min);}
    else
    {
      //SET UPPER BOUND
      double eta = (u2-u1)/varphi2;
      double c2max = 0;
      if (2*M*u2 < 1)
        {c2max += eta*eta;}
      else
        {c2max += (u2+eta)*(u2+eta);}
      // check upper bound is on other side of solution. Since deltaphi_max
      // is positive, this has to be negative. If not, switch limits and check
      // again.
      while (BH_int_dvarphi_du(u1, u2, M, c2max) > 0)
      {
        c2min = c2max*1;
        c2max *= 2;
      }
      // SOLVE WITH BISECTION
      auto BH_tosolve = [u1, u2, varphi2, M](double c2)
                      {return BH_int_dvarphi_du(u1, u2, M, c2)-varphi2;};
      return std::sqrt(bisection(BH_tosolve, c2min, c2max, 1e-5));
    }
  }
  // U2<U1 -> USE MINUS (in dvarphi_du it means switch u1 and u2)
  else 
  {
    double deltaphi_max = BH_int_dvarphi_du(u2, u1, M, c2min);
    if (deltaphi_max <= 0)
      {return std::sqrt(c2min);}
    else
    {
      //SET UPPER BOUND
      double eta = (u1-u2)/varphi2;
      double c2max = 0;
      if (2*M*u1 < 1)
        {c2max += eta*eta;}
      else
        {c2max += (u1+eta)*(u1+eta);}
      // check upper bound is on other side of solution. Since deltaphi_max
      // is positive, this has to be negative. If not, switch limits and check
      // again.
      while (BH_int_dvarphi_du(u2, u1, M, c2max) > 0)
      {
        c2min = c2max*1;
        c2max *= 2;
      }
      // SOLVE WITH BISECTION
      auto BH_tosolve = [u2, u1, varphi2, M](double c2)
                      {return BH_int_dvarphi_du(u2, u1, M, c2)-varphi2;};
      return std::sqrt(bisection(BH_tosolve, c2min, c2max, 1e-5));
    }
  }
};
double BH_time_caus_check(double u1, double u2, double t2, double M, double c)
  {return BH_int_dt_du (u1, u2, M, c) <= t2;}


int main()
{

    // cout<<"\n============== TESTING BISECTION ===============\n";
    // auto r = boost::math::tools::bisect(poly, 0, 8, TermCond());
    // cout<<"Boost's solution to x^3-8 is "<<(r.first+r.second)/2<<endl;

    cout<<"\n================= TESTING BH DIFF EQS ==================\n";
    cout<<"\n====== 1. Testing the integral of d(varphi)/du ======\n";
    double M = 1;

    cout<<"Test 1.1\n";
    double c = 0.5;
    double u1 = 1.;
    double u2 = 2.;
    cout<<"M=1; c=0.5; u1=1; u2=2        -> W.Alfa: 0.5002\n";
    cout<<"                              -> We say: "<<BH_int_dvarphi_du
                                                      (u1,u2,M,c*c)
        <<endl;
    
    cout<<"Test 1.2\n";
    c  = 1.25;
    u1 = 1./3.;
    u2 = 2.;
    cout<<"M=1; c=1.25; u1=1/3; u2=2     -> W.Alfa: 0.9194\n";
    cout<<"                              -> We say: "<<BH_int_dvarphi_du
                                                       (u1,u2,M,c*c)
        <<endl;

    cout<<"Test 1.3\n";
    c = 0.5;
    u1 = 1./2.;
    u2 = 1./3.;
    cout<<"M=1; c=0.5; u1=1/2; u2=1/3    -> W.Alfa: 0.6522\n";
    cout<<"                              -> We say: "<<BH_int_dvarphi_du
                                                        (u1,u2,M,c*c)
        <<endl;
    

    cout<<"\n========= 2. Testing the solver for c^2 =========\n";

    cout<<"Test 1.1\n";
    double varphi2 = 0.5002;
    u1 = 1.;
    u2 = 2.;
    cout<<"M=1; vphi=0.5002; u1=1; u2=2  -> W.Alfa: 0.25\n";
    cout<<"                              -> We say: "<<BH_c2_solver
                                                        (u1, u2, varphi2, M)
        <<endl;
    
    cout<<"Test 1.2\n";
    varphi2 = 0.9194;
    u1 = 1./3.;
    u2 = 2.;
    cout<<"M=1; vphi=0.5002; u1=1/3; u2=2-> W.Alfa: 1.5625\n";
    cout<<"                              -> We say: "<<BH_c2_solver
                                                       (u1, u2, varphi2, M)
        <<endl;

    cout<<"Test 1.3\n";
    varphi2 = 0.6522;
    u1 = 1./2.;
    u2 = 1./3.;
    cout<<"M=1; varphi2 = 0.6522; u1=1/2; u2=1/3-> W.Alfa: 0.25\n";
    cout<<"                                     -> We say: "<<BH_c2_solver
                                                           (u1, u2, varphi2, M)
        <<endl;


    cout<<"\n======= 3. Testing the integral of dt/du =======\n";

    cout<<"Test 3.1\n";
    c = 0.75;
    u1 = 1./3;
    u2 = 2.;
    cout<<"M=1; c=0.75; u1=1/3; u2=2     -> W.Alfa: 3.15412\n";
    cout<<"                              -> We say: "<<BH_int_dt_du
                                                       (u1,u2,M,c*c)
        <<endl;
    
    cout<<"Test 3.2\n";
    c  = 1.25;
    u1 = 1./3.;
    u2 = 1/2.;
    cout<<"M=1; c=1.25; u1=1/3; u2=1/2   -> W.Alfa: 1.054\n";
    cout<<"                              -> We say: "<<BH_int_dt_du
                                                       (u1,u2,M,c*c)
        <<endl;

    cout<<"Test 3.3\n";
    c = 0.5;
    u1 = 3;
    u2 = 1.25;
    cout<<"M=1; c=0.5; u1=3; u2=1.25     -> W.Alfa: -0.695754\n";
    cout<<"                              -> We say: "<<BH_int_dt_du
                                                        (u1,u2,M,c*c)
        <<endl;
    
    cout<<"Test 3.4\n";
    c = 1.25;
    u1 = 0.666;
    u2 = 0.5;
    cout<<"M=1; c=1.25; u1=0.666; u2=0.5 -> W.Alfa: 0.6522\n";
    cout<<"                              -> We say: "<<BH_int_dvarphi_du
                                                        (u1,u2,M,c*c)
        <<endl;
    

    return 0;
}