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
#include <map>
//#include <boost/asio.hpp>

#include "functions.h"
#include "vecfunctions.h"
//#include "causet.h"
//#include "embeddedcauset.h"
//#include "shapes.h"
//#include "spacetimes.h"


using namespace std::chrono;
using namespace std;
int main(){


    ////TESTING MAPS' METHODS
    /////////////////////////
        // std::map<const char*, double> mymap = {{"a", 1},{"b",2}, {"c",3}};

        // mymap["c"] = 10;
        // std::cout<< "After 'mymap['c'] = 10', the size is "<<mymap.size()<<endl;
        // double a = mymap["c"];
        // std::cout<< "And the element is "<<a<<endl;
        // std::cout<< "After assignment, the size is "<<mymap.size()<<endl;
        // a++;
        // std::cout<< "After modification of a, the element is "<<mymap["c"]<<endl;

        // mymap.insert({"c", 100});
        // std::cout<< "After 'mymap.insert({'c', 100})', the size is "<<mymap.size()<<endl;
        // double b = mymap.find("c")->second;
        // std::cout<< "and the element is "<<b<<endl;
        // std::cout<< "After assignment, the size is "<<mymap.size()<<endl;
        // b++;
        // std::cout<< "After modification of b, the element is "<<mymap["c"]<<endl;
        // mymap["d"] = 5;
        // std::cout<< "After 'mymap['d'] = 5', the size is "<<mymap.size()<<endl;
    

    // TESTING INTERSECTION OF SETS
    ////////////////////////////////
        unordered_set<int> a = {1,2,3,4,5,6,7,8,9,10};
        unordered_set<int> b = {1,3,5,6,8,10,14,22,37,48,59,62};
        cout<<"Sets to intersect:\n";
        cout<<"a = { ";
        for (auto a_i : a)
            {cout<<a_i<<", ";}
        cout<<"}\n";
        cout<<"b = { ";
        for (auto b_i : b)
            {cout<<b_i<<", ";}
        cout<<"}\n";
        cout<<"Intersection:\n";
        auto result = set_intersection(a, b);
        cout<<"r = {";
        for (auto x : result)
            {cout<<x<<", ";}
        cout<<"}\n";


// int count = 3;
// int dim = 2;
// vector<vector<double>> coords(count,(vector<double>(dim,0.0)));
// std::cout << coords.size() << std::endl;

// for (int i=0; i<count; i++)
// {
//     for (int j=0; j<count; j++)
//     {
//         coords[i][j] = 1.0;
//     }
// }
// print_vector(coords);

// std::cout << "Hello World" << std::endl;

// std::set<int> s1 = {1,2,3,4};
// std::set<int> s2 = {2,4};

auto start = high_resolution_clock::now();

// int k = 1;
// for (int i=1; i<1000000001; i++)
// {
//     k=i;
//     //std::cout << k << std::endl;
// }

// std::cout << "k =" << k << std::endl;

// std::set<int> s = set_intersection(s1,s2);
// print_set(s);

//// TEST RANDOM INTERGERS
///////////////////////////
// std::vector<int> a;
// a = distinct_randint(100000,1000000);
// //print_vector(a);
// std::cout << a[1000] << std::endl;
/*
int DIM = 4;
int N = 10000;
vector<vector<double>> coords = generate_2Dvector(N,DIM,0,2);
//std::cout << "This file works" << std::endl;



Causet c(coords,"cmatrix");
*/
auto stop = high_resolution_clock::now();
double duration = duration_cast<microseconds>(stop - start).count();
std::cout << "Time taken: "
         << duration/pow(10,6) << " seconds" << std::endl;


};  
