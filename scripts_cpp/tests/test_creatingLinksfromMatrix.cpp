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


#include "../causets_cpp/functions.h"
#include "../causets_cpp/vecfunctions.h"

#include <boost/range/combine.hpp>
#include <omp.h>

using namespace std::chrono;



// Variables to change
int _size = 10000;
double min = 0;
double max = 1;




int main(){


// Declerations
std::vector<std::unordered_set<int>> _future_links;
_future_links.resize(_size);
std::vector<std::vector<int>> _CMatrix = generate_2Dvector_of_ints(_size, _size,
                                                                    min, max);

auto beginning = high_resolution_clock::now();

#pragma omp parallel for
for (int i=0; i<_size; i++)
{
    for (int j=i+1; j<_size; j++)
    {
        if (_CMatrix[i][j] == 0) {
            continue;
        }
        else
        {
            bool has_broken = false;
            for (int k=i+1; k<j;k++)
            {
                if (_CMatrix[i][k]*_CMatrix[k][j]!=0){
                    has_broken = true;
                    break;}
            }
            if (!has_broken){
                _future_links[i].insert(j);}
        }
    }
}

auto finish = high_resolution_clock::now();
double duration = duration_cast<microseconds>(finish - beginning).count();
std::cout << "\nProgram took in total: "
        << duration/pow(10,6) << " seconds\n" << std::endl;


}


