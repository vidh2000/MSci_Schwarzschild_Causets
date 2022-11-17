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

#include "causets_cpp/embeddedcauset.h"
#include "causets_cpp/shapes.h"
#include "causets_cpp/spacetimes.h"
#include "causets_cpp/functions.h"
#include "causets_cpp/vecfunctions.h"

#include "boost/math/tools/roots.hpp"
#include "boost/numeric/odeint.hpp"

using namespace boost::numeric::odeint;
using std::cout;
using std::endl;
using std::vector;

int main()
{
    cout<<"\n================= TESTING BH DIFF EQS ==================\n";
    double c, u1, u2, varphi2;
    double M = 1;

    cout<<"\n\n====== 1. Testing the integral of d(varphi)/du ======\n";
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
    

    cout<<"\n\n========= 2. Testing the solver for c =========\n";

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
    varphi2 = 0.352;
    u1 = 1./2.;
    u2 = 1./3.;
    cout<<"M=1; varphi2 = 0.352; u1=1/2; u2=1/3-> W.Alfa: 0.5\n";
    cout<<"                                     -> We say: "
        <<Spacetime::BH_c_solver(u1, u2, varphi2, M)<<endl;


    cout<<"\n\n======= 3. Testing the integral of dt/du =======\n";

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
    

    cout<<"\n\n============= TESTING POINTS CAUSALITY IN BH ==============\n";
    cout<<"Use ex in 'He, Rideout. Causal Sets Black Hole. 2009. page 10'.\n";

    //Prepare 4D Spacetime(BlackHole) & Shape(whatever)
    Spacetime St = Spacetime();
    St.BlackHoleSpacetime(4, 1);
    CoordinateShape Sh = CoordinateShape(4);


    //Take  EF coords from "He, Rideout. Causal Sets Black Hole. Page 10."
    vector<vector<double>> coords = {
                                {0.410895, 2.36161,  1.80295, 0.57951},//out
                                {1.109415, 2.89891,  1.04335, 4.25531},//out
                                {1.133105, 1.36083,  1.89919, 1.06482},//in
                                {2.743428, 2.74093,  2.97906, 4.22204},//out
                                {3.235970, 0.65462,  0.11664, 5.06884},//in
                                {3.972871, 0.96354,  2.33727, 1.38169},//in
                                {5.230757, 2.34476,  1.11855, 3.47242},//out
                                {6.014261, 0.664739, 2.82235, 0.95459}, //in
                                {6.193089, 0.429636, 2.20122, 1.99644}}; //in
    
    //Define Embedded Causet (just the causet, no causal relations)
    EmbeddedCauset C = EmbeddedCauset(St, Sh, coords, false, false, false,
                                      false, false, "past");
    

    //In depth analysis of all generic
    cout<<"\n\n======== TEST GENERIC CONNECTIONS =========\n";
    // double u1, u2, t1, t2, varphi2;
    // double M = 1;
    cout<<"\nChecking 0 -> 4: (should have c^2=0.04604, geo_time=4.6979) \n";
    bool c04 = C.AprecB(coords[0], coords[4]);
    cout<<"Should give 0 -> "<<c04<<endl;
    cout<<"\nChecking 0 -> 5: (should have c^2=0.2727, geo_time=2.119) \n";
    bool c05 = C.AprecB(coords[0], coords[5]); 
    cout<<"Should give 1 -> "<<c05<<endl;
    cout<<"\nChecking 1 -> 5: (should have c^2=0.0388, geo_time=9.8635) \n";
    bool c15 = C.AprecB(coords[1], coords[5]); 
    cout<<"Should give 0 -> "<<c15<<endl;
    cout<<"\nChecking 2 -> 5: (should get to the integral and get c=-1)\n";
    bool c25 = C.AprecB(coords[2], coords[5]);
    cout<<"Should give 0 -> "<<c25<<endl;
    cout<<"\nChecking 1 -> 6: (should have c^2=0.0476, geo_time=2.6097)\n";
    bool c16 = C.AprecB(coords[1], coords[6]);
    cout<<"Should give 1 -> "<<c16<<endl;
    cout<<"\nChecking 3 -> 7: (should have c^2=4.55, geo_time=2.1833)\n";
    bool c37 = C.AprecB(coords[3], coords[7]);
    cout<<"Should give 1 -> "<<c37<<endl;
    cout<<"\nChecking 5 -> 8: (should have c^2=1.356, geo_time=0.6671)\n";
    bool c58 = C.AprecB(coords[5], coords[8]);
    cout<<"Should give 1 -> "<<c58<<endl;

    //Compare
    cout<<"\n\n======== TEST CAUSALITIES OF THE 9 POINTS =========\n";
    //Compute Causality Matrix & Sets
    C.make_all_but_links();
    cout<<"Their Causality Matrix is (5 means unspecified): \n";
    vector<vector<double>> HRmatrix = {
                                        {0,0,0, 5,0,1, 0,5,1},
                                        {0,0,0, 5,0,0, 1,5,5},
                                        {0,0,0, 0,0,0, 0,5,5},

                                        {0,0,0, 0,5,5, 5,1,5},
                                        {0,0,0, 0,0,0, 0,5,5},
                                        {0,0,0, 0,0,0, 0,5,1},

                                        {0,0,0, 0,0,0, 0,5,5},
                                        {0,0,0, 0,0,0, 0,0,5},
                                        {0,0,0, 0,0,0, 0,0,0},
                                      };
    print_vector(HRmatrix);
    cout<<"\nOur Causality Matrix is: \n";
    print_vector(C._CMatrix);
    cout<<"The difference is : \n";
    cout<<"Note: 0 -> ok, 1 -> should have been related, \
        \n      -1 ->shouldn't have been related, 4, 5 -> unspecified.\n";
    vector<vector<double>> difference (coords.size(), 
                                       vector<double> (coords.size()));
    for (int i = 0; i<coords.size(); i++)
    {
        for (int j = 0; j<coords.size(); j++)
            {difference[i][j] = HRmatrix[i][j] - C._CMatrix[i][j];}
    }
    print_vector(difference);
    // Also print future sets, cause, why not
    for (int i = 0; i<coords.size(); i++)
    {
        cout<<"\nFuture of "<<i<<": ";
        for (int ej : C._futures[i])
            {cout<<ej<<" ";}
    }
    
    


    return 0;
}