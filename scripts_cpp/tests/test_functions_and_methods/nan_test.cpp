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



using namespace std::chrono;
using namespace std;
int main(){
    auto start = high_resolution_clock::now();

 
    cout<<"\n===============TESTING NANS METHODS===========\n";
    double a = std::nan("");
    double b = 3;

    std::cout<<"a = "<<a<<std::endl;
    std::cout<<"b>a?"<<(b> a)<<std::endl;
    std::cout<<"b<a?"<<(b< a)<<std::endl;
    std::cout<<"b=a?"<<(b==a)<<std::endl;
    

    auto stop = high_resolution_clock::now();
    double duration = duration_cast<microseconds>(stop - start).count();
    std::cout << "Time taken: "
            << duration/pow(10,6) << " seconds" << std::endl;


};  