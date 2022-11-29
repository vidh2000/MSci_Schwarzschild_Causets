#include <algorithm>
#include <chrono>
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
#include <unordered_set>

using namespace std::chrono;
using std::vector;
using std::cout;
using std::endl;
using std::unordered_set;


int main(){

    cout<<"\n============= PUSH_BACK VS INSERT ==================\n";
    int imax = 1e5;

    auto start = high_resolution_clock::now();
    for (int i = 0; i<imax; i++){    
    }
    auto end = high_resolution_clock::now();
    double duration = duration_cast<microseconds>(end - start).count();
    cout << "EMPTY Looping   = " << duration/pow(10,6)<< " s" << endl;

    start = high_resolution_clock::now();
    vector<int> a;
    for (int i = 0; i<1e5; i++){
        a.push_back(i);
    }
    auto vend = high_resolution_clock::now();
    duration = duration_cast<microseconds>(vend - start).count();
    cout << "VECTOR Push Back = " << duration/pow(10,6)<< " s" << endl;

    start = high_resolution_clock::now();
    unordered_set<int> b;
    for (int i = 0; i<1e5; i++){
        b.insert(i);
    }
    auto usend = high_resolution_clock::now();
    duration = duration_cast<microseconds>(usend - start).count();
    cout << "UNORD_SET Insert = " << duration/pow(10,6)<< " s" << endl;


    cout<<"\n============= VECTOR VS ARRAY ==================\n";
    int size = 1e5; 

    auto VA_start = high_resolution_clock::now();
    vector<vector<int>> V;
    V.resize(size, vector<int>(size,0));
    for (int i = 0; i<V.size(); i++){
        for (int j = i+1; j<V.size(); j++){
            V[i][j] = 1;}
    }
    auto Vend = high_resolution_clock::now();
    duration = duration_cast<microseconds>(Vend - VA_start).count();
    cout << "VECTOR Matrix Looping = " << duration/pow(10,6)<< " s" << endl;

    VA_start = high_resolution_clock::now();
    int A[size][size];
    for (int i = 0; i<sizeof(A)/sizeof(A[0]); i++){
        for (int j = i+1; j<sizeof(A[i])/sizeof(A[i][0]); j++){
            A[i][j] = 1;}
    }
    auto Aend = high_resolution_clock::now();
    duration = duration_cast<microseconds>(Aend - VA_start).count();
    cout << "ARRAY Matrix Looping  = " << duration/pow(10,6)<< " s" << endl;

    return 0;

}