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
using namespace std::chrono;

#include "../scripts_cpp/causets_cpp/functions.h"
#include "../scripts_cpp/causets_cpp/vecfunctions.h"

#include <omp.h>



int repeat = 1000;
int n = 1000;


int N_cores;


int main(){

// enable/disable OpenMP via include in tasks.json
#pragma omp parallel
N_cores = omp_get_thread_num()+1;


std::cout << "Number of CPU cores= " << N_cores << std::endl;


// No parallelism
std::vector<std::vector<double>> a = generate_2Dvector(n,n,0,n);
a = generate_2Dvector(n,n,0,n);

auto start = high_resolution_clock::now();

for (int r=0; r<repeat; r++)
{
    for(int i=0; i<n; i++) {
        for(int j=0; j<n; j++) {
            a[i][j] += j+i;
        }
    }
}
auto stop = high_resolution_clock::now();
auto duration = duration_cast<microseconds>(stop - start).count();
std::cout << "Without parallelism: "
        << duration/pow(10,6) << " seconds" << std::endl;


// Method 1
a = generate_2Dvector(n,n,0,n);
auto start1 = high_resolution_clock::now();


for (int r=0; r<repeat; r++)
{
    #pragma omp parallel for schedule(dynamic,N_cores)
    for(int i=0; i<n; i++) {
        for(int j=0; j<n; j++) {
            a[i][j] += j+i;
        }
    }
}
auto stop1 = high_resolution_clock::now();
auto duration1 = duration_cast<microseconds>(stop1 - start1).count();
std::cout << "Method 1: "
        << duration1/pow(10,6) << " seconds" << std::endl;

// Method 2
a = generate_2Dvector(n,n,0,n);
auto start2 = high_resolution_clock::now();

for (int r=0; r<repeat; r++)
{
    #pragma omp for nowait
    for(int i=0; i<n; i++) {
        for(int j=0; j<n; j++) {
            a[i][j] += j+i;
        }
    }
}
auto stop2 = high_resolution_clock::now();
auto duration2 = duration_cast<microseconds>(stop2 - start2).count();
std::cout << "Method 2: "
        << duration2/pow(10,6) << " seconds" << std::endl;

// Method 3
a = generate_2Dvector(n,n,0,n);
auto start3 = high_resolution_clock::now();

for (int r=0; r<repeat; r++)
{
    #pragma omp parallel for schedule(dynamic,N_cores)
    for(int i=0; i<n; i++) {
        for(int j=0; j<n; j++) {
            a[i][j] += j+i;
        }
    }
}
auto stop3 = high_resolution_clock::now();
auto duration3 = duration_cast<microseconds>(stop3 - start3).count();
std::cout << "Method 3: "
        << duration3/pow(10,6) << " seconds" << std::endl;

// Method 4
a = generate_2Dvector(n,n,0,n);
auto start4 = high_resolution_clock::now();


for (int r=0; r<repeat; r++)
{
    
    for(int i=0; i<n; i++) {
        for(int j=0; j<n; j++) {
            a[i][j] += j+i;
        }
    }
}
auto stop4 = high_resolution_clock::now();
auto duration4 = duration_cast<microseconds>(stop4 - start4).count();
std::cout << "Method 4: "
        << duration4/pow(10,6) << " seconds" << std::endl;


return 0;
}



