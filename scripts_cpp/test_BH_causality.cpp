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
#include "causets_cpp/spacetimes.h"
#include "boost/math/tools/roots.hpp"
#include "boost/numeric/odeint.hpp"

using namespace boost::numeric::odeint;
using std::cout;
using std::endl;

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
    cout<<"                              -> We say: "
        <<Spacetime::BH_int_dvarphi_du(u1,u2,c*c,M)<<endl;
    
    cout<<"Test 1.2\n";
    c  = 1.25;
    u1 = 1./3.;
    u2 = 2.;
    cout<<"M=1; c=1.25; u1=1/3; u2=2     -> W.Alfa: 0.919415\n";
    cout<<"                              -> We say: "
        <<Spacetime::BH_int_dvarphi_du(u1,u2,c*c,M)<<endl;

    cout<<"Test 1.3\n";
    c = 0.5;
    u1 = 1./2.;
    u2 = 1./3.;
    cout<<"M=1; c=0.5; u1=1/2; u2=1/3    -> W.Alfa: 0.352028\n";
    cout<<"                              -> We say: "
        <<Spacetime::BH_int_dvarphi_du(u1,u2,c*c,M)<<endl;
    

    cout<<"\n========= 2. Testing the solver for c =========\n";

    cout<<"Test 2.1\n";
    varphi2 = 0.500168;
    u1 = 1.;
    u2 = 2.;
    cout<<"M=1; vphi=0.5002; u1=1; u2=2  -> W.Alfa: 0.5\n";
    cout<<"                              -> We say: "
        <<Spacetime::BH_c_solver(u1, u2, varphi2, M)<<endl;
    
    cout<<"Test 2.2\n";
    varphi2 = 0.919415;
    u1 = 1./3.;
    u2 = 2.;
    cout<<"M=1; vphi=0.5002; u1=1/3; u2=2-> W.Alfa: 1.25\n";
    cout<<"                              -> We say: "
        <<Spacetime::BH_c_solver(u1, u2, varphi2, M)<<endl;

    cout<<"Test 2.3\n";
    varphi2 = 0.352028;
    u1 = 1./2.;
    u2 = 1./3.;
    cout<<"M=1; varphi2 = 0.6522; u1=1/2; u2=1/3-> W.Alfa: 0.5\n";
    cout<<"                                     -> We say: "
        <<Spacetime::BH_c_solver(u1, u2, varphi2, M)<<endl;


    cout<<"\n======= 3. Testing the integral of dt/du =======\n";

    cout<<"Test 3.1\n";
    c = 0.75;
    u1 = 1./3.;
    u2 = 1.5;
    cout<<"M=1; c=0.75; u1=1/3; u2=1.5   -> W.Alfa: 2.93751\n";
    cout<<"                              -> We say: "
        <<Spacetime::BH_int_dt_du(u1,u2,c,M)<<endl;
    
    cout<<"Test 3.2\n";
    c = 0.75;
    u1 = 0.75;
    u2 = 2;
    cout<<"M=1; c=0.75; u1=0.75; u2=2    -> W.Alfa: 1.13855\n";
    cout<<"                              -> We say: "
        <<Spacetime::BH_int_dt_du(u1,u2,c,M)<<endl;
    
    cout<<"Test 3.3\n";
    c  = 1.25;
    u1 = 1./4.;
    u2 = 1/3.;
    cout<<"M=1; c=1.25; u1=1/4; u2=1/3   -> W.Alfa: 1.02712\n";
    cout<<"                              -> We say: "
        <<Spacetime::BH_int_dt_du(u1,u2,c,M)<<endl;

    cout<<"Test 3.4\n";
    c = 0.5;
    u1 = 3;
    u2 = 1.25;
    cout<<"M=1; c=0.5; u1=3; u2=1.25     -> W.Alfa: -0.695754\n";
    cout<<"                              -> We say: "
        <<Spacetime::BH_int_dt_du(u1,u2,c,M)<<endl;
    
    cout<<"Test 3.5\n";
    c = 1.25;
    u1 = 0.666;
    u2 = 0.55;
    cout<<"M=1; c=1.25; u1=0.666; u2=0.55-> W.Alfa: -3.68198\n";
    cout<<"                              -> We say: "
        <<Spacetime::BH_int_dt_du(u1,u2,c,M)<<endl;
    

    cout<<"\n============= TESTING POINTS CAUSALITY IN BH ==============\n";
    cout<<"Using Example in He, Rideout 2009\n";
    

    return 0;
}