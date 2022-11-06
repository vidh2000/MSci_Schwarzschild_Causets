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

// Test C++ version... comment Vid New
int main() {

std::cout << "Hello World" << std::endl;

std::random_device rd;
int seed = rd();
std::mt19937 gen(seed);
std::uniform_real_distribution<> dis(0,5);
int e1 = (int) dis(gen), e2 =(int) dis(gen);

std::cout << e1 << " " << e2 << std::endl;

/*
int DIM = 4;
int N = 10000;
vector<vector<double>> coords = generate_2Dvector(N,DIM,0,2);
//std::cout << "This file works" << std::endl;

auto start = high_resolution_clock::now();

Causet c(coords,"cmatrix");

auto stop = high_resolution_clock::now();
double duration = duration_cast<microseconds>(stop - start).count();
std::cout << "Time taken by function in D=" << DIM << ": "
         << duration/pow(10,6) << " seconds" << std::endl;
*/

};  
