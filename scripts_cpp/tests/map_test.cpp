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

#include "../causets_cpp/functions.h"
#include "../causets_cpp/vecfunctions.h"
//#include "causet.h"
//#include "embeddedcauset.h"
//#include "shapes.h"
//#include "../causets_cpp/spacetimes.h"

template<typename num>
void print_num_map(std::map<int, num> const &map)
{
    for (auto const &pair: map) {
        std::cout << "{" << pair.first << ": " << pair.second << "}\n";
    }
}


using namespace std::chrono;
using namespace std;
int main(){
    auto start = high_resolution_clock::now();

 
    cout<<"\n===============TESTING MAPS METHODS===========\n";

    std::map<int, int> mymap;
    int a = mymap[1];
    cout<<"\nAfter 'int a = mymap[1]', a is "<<a<<endl;
    cout<<"and the map is"<<endl;
    print_num_map(mymap);
    mymap[1] += 1;
    cout<<"\nAfter 'mymap[1] += 1', mymap is"<<endl;
    print_num_map(mymap);
    mymap[1] += 1;
    cout<<"\nAfter another 'mymap[1] += 1'"<<endl;
    print_num_map(mymap);
    mymap[2] = 0;
    cout<<"\nAfter 'mymap[2] = 0' we have mymap[1] = "<<endl;
    print_num_map(mymap);
    mymap[2] += 1;
    cout<<"\nAfter 'mymap[2] += 1' we have mymap[1] = "<<endl;
    print_num_map(mymap);

    auto stop = high_resolution_clock::now();
    double duration = duration_cast<microseconds>(stop - start).count();
    std::cout << "Time taken: "
            << duration/pow(10,6) << " seconds" << std::endl;


};  