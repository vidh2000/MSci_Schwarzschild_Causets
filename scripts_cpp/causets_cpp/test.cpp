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

std::map<const char*, double> p = {{"ha", 1}};
std::cout << p["ha"] << std::endl;
std::cout << p["nah"] << std::endl;
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


std::vector<int> a;
a = distinct_randint(100000,1000000);
//print_vector(a);
std::cout << a[1000] << std::endl;
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
