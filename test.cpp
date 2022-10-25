#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <fstream>
#include <stack>
#include <string>
#include <stdio.h>
#include <vector>
#include <set>



std::set<int> s;

int main(){

s = {1,2,3,42};

std::cout << "The elements in set are: ";
    for (auto it = s.begin(); it != s.end(); it++)
        std::cout << *it << " ";


}