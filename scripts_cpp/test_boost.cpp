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

#include "causets_cpp/functions.h"
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


// VERY IMPORTANT NOTE: CANNOT COPY-PASTE THESE FUNCTIONS IN SPACETIME
// BECAUSE IN SPACETIME ORDER OF SOME ARGUMENTS IS REVERSED: MASS ALWAYS LAST
// BRINGING CHANGES IN DECLARATION AND CALLS INSIDE OTHER FUNCTIONS
void BH_dvarphi_du (double& dpdu, double u, double M, double c2)
{
  dpdu = std::pow(2*M*u*u*u - u*u + c2, -0.5);
};
double BH_int_dvarphi_du(double u1, double u2, double M, double c2)
{
    auto BH_dvarphi_du_forint = [M, c2]
                            (const double& varphi, double& dpdu, const double u)
                            {dpdu = std::pow(u*u*(2*M*u - 1) + c2, -0.5);};
    double varphi = 0;
    if (u2>=u1)
    {
        boost::numeric::odeint::integrate(BH_dvarphi_du_forint, varphi, 
                                          u1, u2, (u2-u1)/20.);
    }
    else
    {
        boost::numeric::odeint::integrate(BH_dvarphi_du_forint, varphi, 
                                          u2, u1, (u1-u2)/20.);
    }
    return varphi;
};

double BH_c_solver (double u1, double u2, double varphi2, double M)
{
  // U2>=U1 -> USE PLUS
  if (u2>=u1) 
  {
    //SET MINIMUM BOUND
    double c2min = 0;
    if (2*M*u1<=1.)
    {
        double D1 = u1*u1 *(1 - 2*M*u1);
        double D2 = u2*u2 *(1 - 2*M*u2);
        c2min += (D1>D2)? D1 : D2;
        c2min += 1e-3; //for anti-divergence purposes
    }
    //SET UPPER BOUND
    double deltaphi_max = BH_int_dvarphi_du(u1, u2, M, c2min);
    if (deltaphi_max <= varphi2)
      {return std::sqrt(c2min);}
    else
    {
      double eta = (u2-u1)/varphi2;
      double c2max = 0;
      if (2*M*u2 < 1)
        {c2max += eta*eta;}
      else
        {c2max += (u2+eta)*(u2+eta);}
      // check upper bound is on other side of solution. Since deltaphi_max
      // is positive, this has to be negative. If not, double upper limit,
      // turn lower limit to previous upper and check again.
      while (BH_int_dvarphi_du(u1, u2, M, c2max) > varphi2)
      {
        std::cout<<"No biggy, just note in c2solver had to update upper limit"
                 <<std::endl;
        c2min = c2max*1;
        c2max *= 2;
      }
      // SOLVE WITH BISECTION
      auto BH_tosolve = [u1, u2, varphi2, M](double c2)
                      {return BH_int_dvarphi_du(u1, u2, M, c2)-varphi2;};
      return std::sqrt(bisection(BH_tosolve, c2min, c2max, 1e-8));
    }
  }
  // U2<U1 -> USE MINUS (in dvarphi_du it means switch u1 and u2)
  else 
  {
    //SET MINIMUM BOUND
    double c2min = 0;
    if (2*M*u2<=1.)
    {
        double D1 = u1*u1 *(1 - 2*M*u1);
        double D2 = u2*u2 *(1 - 2*M*u2);
        c2min += (D1>D2)? D1 : D2;
        c2min += 1e-3; //for anti-divergence purposes
    }
    //SET UPPER BOUND
    double deltaphi_max = BH_int_dvarphi_du(u2, u1, M, c2min);
    if (deltaphi_max <= 0)
      {return std::sqrt(c2min);}
    else
    {
      double eta = (u1-u2)/varphi2;
      double c2max = 0;
      if (2*M*u1 < 1)
        {c2max += eta*eta;}
      else
        {c2max += (u1+eta)*(u1+eta);}
      // check upper bound is on other side of solution. Since deltaphi_max
      // is positive, this has to be negative. If not, switch limits and check
      // again.
      while (BH_int_dvarphi_du(u2, u1, M, c2max) > varphi2)
      {
        std::cout<<"No biggy, just note in c2solver had to update upper limit"
                 <<std::endl;
        c2min = c2max*1;
        c2max *= 2;
      }
      // SOLVE WITH BISECTION
      auto BH_tosolve = [u2, u1, varphi2, M](double c2)
                      {return BH_int_dvarphi_du(u2, u1, M, c2)-varphi2;};
      return std::sqrt(bisection(BH_tosolve, c2min, c2max, 1e-3));
    }
  }
};
void BH_dt_du_plus (double&dtdu, double u, double M, double c)
{
  double D = u*u*(1-2*M*u);
  dtdu = (c*std::pow(c*c - D, -0.5) - 2*M*u)/D;
};
void BH_dt_du_minus (double&dtdu, double u, double M, double c)
{
  double D = u*u*(1-2*M*u);
  dtdu = (-c*std::pow(c*c - D, -0.5) - 2*M*u)/D;
};
double BH_int_dt_du (double u1, double u2, double M, double c)
{
    double t = 0;
    if (u2>=u1) 
    {
        auto BH_dt_du_forint_plus = [M, c]
                                (const double& t, double& dtdu, const double u)
                                {
                                  if (u==0.5)
                                    {dtdu = 0.488889;}
                                  else
                                  {
                                    double D = u*u*(1-2*M*u);
                                    dtdu = (c/std::sqrt(c*c - D) - 2*M*u)
                                          /D;
                                  }
                                };
        boost::numeric::odeint::integrate(BH_dt_du_forint_plus, t, 
                                              u1, u2, (u2-u1)/20.);
    }
    else /* u1>u2 */
    {
        auto BH_dt_du_forint_minus = [M, c]
                                (const double& t, double& dtdu, const double u)
                                {
                                  double D = u*u*(1-2*M*u);
                                  dtdu = (-c/std::sqrt(c*c - D) - 2*M*u)
                                        /D;
                                };
        if (0.5*M>=u1 || u2>=0.5*M) //0.5*M NOT in [u2, u1]
        {
            //Avoid divergence at 1-2*M*u=0
            u2 += (u2==0.5*M)? 1e-3*M : 0;
            u1 -= (u2==0.5*M)? 1e-3*M : 0;
            //Compute 
            boost::numeric::odeint::integrate(BH_dt_du_forint_minus, t, 
                                              u1, u2, -(u1-u2)/20.);
        }
        else /*0.5*M IN [u2, u1]*/
        {
            //Compute in 2 steps to avoid divergence
            boost::numeric::odeint::integrate(BH_dt_du_forint_minus, t, 
                                            u1, (0.5+0.0001)*M, -(u1-0.5)/20.);
            boost::numeric::odeint::integrate(BH_dt_du_forint_minus, t, 
                                            (0.5-0.0001)*M, u2, -(0.5-u2)/20.);
        }
    }
    return t;
};
double BH_time_caus_check(double u1, double u2, double t2, double M, double c)
  {return BH_int_dt_du (u1, u2, M, c) <= t2;}


int main()
{
    cout<<"\n================= TESTING BH DIFF EQS ==================\n";
    double c, u1, u2, varphi2;
    double M = 1;

    cout<<"\n====== 1. Testing the integral of d(varphi)/du ======\n";
    M = 1;

    cout<<"Test 1.1\n";
    c = 0.5;
    u1 = 1.;
    u2 = 2.;
    cout<<"M=1; c=0.5; u1=1; u2=2        -> W.Alfa: 0.500168\n";
    cout<<"                              -> We say: "<<BH_int_dvarphi_du
                                                      (u1,u2,M,c*c)
        <<endl;
    
    cout<<"Test 1.2\n";
    c  = 1.25;
    u1 = 1./3.;
    u2 = 2.;
    cout<<"M=1; c=1.25; u1=1/3; u2=2     -> W.Alfa: 0.919415\n";
    cout<<"                              -> We say: "<<BH_int_dvarphi_du
                                                       (u1,u2,M,c*c)
        <<endl;

    cout<<"Test 1.3\n";
    c = 0.5;
    u1 = 1./2.;
    u2 = 1./3.;
    cout<<"M=1; c=0.5; u1=1/2; u2=1/3    -> W.Alfa: 0.352028\n";
    cout<<"                              -> We say: "<<BH_int_dvarphi_du
                                                        (u1,u2,M,c*c)
        <<endl;
    

    cout<<"\n========= 2. Testing the solver for c =========\n";

    cout<<"Test 2.1\n";
    varphi2 = 0.500168;
    u1 = 1.;
    u2 = 2.;
    cout<<"M=1; vphi=0.5002; u1=1; u2=2  -> W.Alfa: 0.5\n";
    cout<<"                              -> We say: "<<BH_c_solver
                                                        (u1, u2, varphi2, M)
        <<endl;
    
    cout<<"Test 2.2\n";
    varphi2 = 0.919415;
    u1 = 1./3.;
    u2 = 2.;
    cout<<"M=1; vphi=0.5002; u1=1/3; u2=2-> W.Alfa: 1.25\n";
    cout<<"                              -> We say: "<<BH_c_solver
                                                       (u1, u2, varphi2, M)
        <<endl;

    cout<<"Test 2.3\n";
    varphi2 = 0.352028;
    u1 = 1./2.;
    u2 = 1./3.;
    cout<<"M=1; varphi2 = 0.3522; u1=1/2; u2=1/3-> W.Alfa: 0.5\n";
    cout<<"                                     -> We say: "<<BH_c_solver
                                                           (u1, u2, varphi2, M)
        <<endl;


    cout<<"\n======= 3. Testing the integral of dt/du =======\n";

    cout<<"Test 3.1\n";
    c = 0.75;
    u1 = 1./3.;
    u2 = 1.5;
    cout<<"M=1; c=0.75; u1=1/3; u2=1.5   -> W.Alfa: 2.93751\n";
    cout<<"                              -> We say: "<<BH_int_dt_du
                                                       (u1,u2,M,c)
        <<endl;
    
    cout<<"Test 3.2\n";
    c = 0.75;
    u1 = 0.75;
    u2 = 2;
    cout<<"M=1; c=0.75; u1=0.75; u2=2    -> W.Alfa: 1.13855\n";
    cout<<"                              -> We say: "<<BH_int_dt_du
                                                       (u1,u2,M,c)
        <<endl;
    
    cout<<"Test 3.3\n";
    c  = 1.25;
    u1 = 1./4.;
    u2 = 1/3.;
    cout<<"M=1; c=1.25; u1=1/4; u2=1/3   -> W.Alfa: 1.02712\n";
    cout<<"                              -> We say: "<<BH_int_dt_du
                                                       (u1,u2,M,c)
        <<endl;

    cout<<"Test 3.4\n";
    c = 0.5;
    u1 = 3;
    u2 = 1.25;
    cout<<"M=1; c=0.5; u1=3; u2=1.25     -> W.Alfa: -0.695754\n";
    cout<<"                              -> We say: "<<BH_int_dt_du
                                                        (u1,u2,M,c)
        <<endl;
    
    cout<<"Test 3.5\n";
    c = 1.25;
    u1 = 0.666;
    u2 = 0.55;
    cout<<"M=1; c=1.25; u1=0.666; u2=0.55-> W.Alfa: -3.68198\n";
    cout<<"                              -> We say: "<<BH_int_dt_du
                                                        (u1,u2,M,c)
        <<endl;
    

    return 0;
}