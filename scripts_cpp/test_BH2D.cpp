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