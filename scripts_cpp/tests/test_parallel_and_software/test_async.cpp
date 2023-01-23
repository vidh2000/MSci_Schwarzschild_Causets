
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

// async
#include <future>

void myFunc()
{
    std::cout << "So, we're in!" << std::endl;
}

double square(double x)
{
    return x*x;
}


static void findResult(std::vector<double>& numbers){}

int main(){

std::vector<double> arr = {1,2,3,4,5,6};


// One by one
std::cout << "One by one...\n";
for (const auto& n : arr)
{
    double result = square(n); 
    std::cout << result << "\n";
}

// Async
std::cout << "Async...\n";
for (const auto& n: arr)
{

}


// end of main
}