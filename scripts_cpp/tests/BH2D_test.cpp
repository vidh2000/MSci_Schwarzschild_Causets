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

#include "../causets_cpp/embeddedcauset.h"
#include "../causets_cpp/sprinkledcauset.h"
#include "../causets_cpp/shapes.h"
#include "../causets_cpp/spacetimes.h"
#include "../causets_cpp/functions.h"
#include "../causets_cpp/vecfunctions.h"

#include "boost/math/tools/roots.hpp"
#include "boost/numeric/odeint.hpp"

using namespace boost::numeric::odeint;
using std::cout;
using std::endl;
using std::vector;


////////////////////////////////
int card = 50;

int main()
{ 
    cout<<"\n\n============= TESTING POINTS CAUSALITY IN BH ==============\n";
    cout<<"Use ex in 'He, Rideout. Causal Sets Black Hole. 2009. page 10'.\n";

    //Prepare 4D Spacetime(BlackHole) & Shape(whatever)
    Spacetime St = Spacetime();
    St.BlackHoleSpacetime(2, 1, "EFuv");
    CoordinateShape Sh = CoordinateShape(2);
    
    //Define Embedded Causet (just the causet, no causal relations)
    SprinkledCauset C = SprinkledCauset(card, St, Sh, true, false, false,
                                        false, false, "past");
    

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