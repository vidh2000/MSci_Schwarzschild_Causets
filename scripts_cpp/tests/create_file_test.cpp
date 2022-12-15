#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <filesystem>
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


int main()
{
    auto path = std::filesystem::current_path();
    std::cout<<path<<std::endl;

    std::ofstream myout1;
    myout1.open("AAAstupid file1.txt");
    std::cout<<"Is myout1 file open? "<<myout1.is_open()<<std::endl;
    myout1.close();
    std::cout<<"Done"<<std::endl;

    std::ofstream myout2("AAAstupid file2.txt");
    std::cout<<"Is myout2 file open? "<<myout2.is_open()<<std::endl;
    myout2.close();
    std::cout<<"Done"<<std::endl;

    std::fstream myout3("AAAstupid file3.txt");
    std::cout<<"Is myout3 file open? "<<myout3.is_open()<<std::endl;
    myout3.close();
    std::cout<<"Done"<<std::endl;
    return 0;
}