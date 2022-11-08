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
//#include <boost/asio.hpp>

#include "functions.h"
#include "vecfunctions.h"
//#include "causet.h"
//#include "embeddedcauset.h"
//#include "shapes.h"
//#include "spacetimes.h"

using namespace std::chrono;

int main(){

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

int b=0;
{std::cout<<!b<<std::endl;}

std::vector<int> a;
a = distinct_randint(2500,10000);
print_vector(a);

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
