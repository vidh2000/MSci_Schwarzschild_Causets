#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
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

#include <unistd.h>
#include <sys/types.h>
//#include <pwd.h>

#include "../../causets_cpp/sprinkledcauset.h"
#include "../../causets_cpp/shapes.h"
#include "../../causets_cpp/spacetimes.h"

#include "../../causets_cpp/functions.h"
#include "../../causets_cpp/vecfunctions.h"

using namespace std::chrono;
using std::cout;
using std::endl;
using std::vector;

//running parameters
bool run_on_hpc = false;

// Sprinkled causet parameters
int myreps = 20;
std::vector<int> dims = {2,3,4};
std::vector<int> mycards = {128, 256, 512, 1024, 2048, 4096, 8192, 16384};
double radius = 1.0;

bool poisson = true;
bool make_matrix = true;
bool special = false;
bool use_transitivity = false;
bool make_sets = false;
bool make_links = false;
const char* sets_type = "futures"; // "both only", "all", "pasts", "futures"

int main(int argc, char** argv) {

    int Nreps = 0;
    std::vector<int> cards;

    if (argc == 1)
    {
        Nreps += myreps;
        cards = mycards;
    }
    else if (argc == 3)
    {
      // Convert the first argument to an integer
      Nreps += std::atoi(argv[1]);
      //Parse the second argument as vector of integers
      std::vector<int> cards;
      char* endptr = nullptr;
      long int value = std::strtol(argv[2], &endptr, 10);
      if (*endptr != '\0') {
          std::cerr << "Error: argument 2 is not a valid integer vector" 
                    << std::endl;
          return 1;
      }
      cards.push_back(static_cast<int>(value));
    }
    else 
    {
        std::cerr << "Must use extra args: <integer> <vector1>" << std::endl;
        return 1;
    }



    auto start = high_resolution_clock::now();
    
    cout<<"\n==============================================\n";
    cout<<"========FIRST TEST ORDERING FRACTION==========\n";
    cout<<"==============================================\n";

    cout<<"\nTEST1.1: 1D CHAIN CAUSET\n";
    int N = 8;
    vector<vector<int>> Cm (N, vector<int> (8));
    for (int i = 0; i<N-1; i++)
    {
        for (int j = i+1; j<N; j++) Cm[i][j] = 1;
    }
    Causet C = Causet(Cm);
    C.make_sets_fromC();
    cout<<"A1 has ord_fr from CMatrix = "<<C.ord_fr(3,6)<<" = 1?\n";
    cout<<"A1 has ord_fr from sets    = "<<C.ord_fr(3,6,"choose",false)<<" = 1?\n";
    cout<<"A1 has ord_fr from CMatrix = "<<C.ord_fr(0,N-1)<<" = 1?\n";
    cout<<"A1 has ord_fr from sets    = "<<C.ord_fr(3,6,"choose",false)<<" = 1?\n";


    cout<<"\nTEST1.2: A GENERAL CAUSET\n";
    /* 
      IT LOOKS LIKE THIS:
          6   8  9  10
          \  / \ | /
            5    7
             \  /
               4
               |
               3
             /   \
            1     2
             \   /
               0
    */
    Cm = {{0, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1},
          {0, 0, 0,  1, 1, 1,  1, 1, 1,  1, 1},
          {0, 0, 0,  1, 1, 1,  1, 1, 1,  1, 1},

          {0, 0, 0,  0, 1, 1,  1, 1, 1,  1, 1},
          {0, 0, 0,  0, 0, 1,  1, 1, 1,  1, 1},
          {0, 0, 0,  0, 0, 0,  1, 0, 1,  0, 0},

          {0, 0, 0,  0, 0, 0,  0, 0, 0,  0, 0},
          {0, 0, 0,  0, 0, 0,  0, 0, 1,  1, 1},
          {0, 0, 0,  0, 0, 0,  0, 0, 0,  0, 0},

          {0, 0, 0,  0, 0, 0,  0, 0, 0,  0, 0},
          {0, 0, 0,  0, 0, 0,  0, 0, 0,  0, 0}};
    C = Causet(Cm);
    cout<<"Has ord_fr from CMatrix = "<< C.ord_fr(3,8)<<endl;
    cout<<"Has ord_fr from sets    = "<< C.ord_fr(3,8,"choose",true)<<endl;
    cout<<"Should have 0.9\n";
    cout<<"Has ord_fr from CMatrix = "<< C.ord_fr(2,8)<<endl;
    cout<<"Has ord_fr from sets    = "<< C.ord_fr(2,8,"choose",true)<<endl;
    cout<<"Should have 0.93\n";
    cout<<"Has ord_fr from CMatrix = "<< C.ord_fr(3,9)<<endl;
    cout<<"Has ord_fr from sets    = "<< C.ord_fr(3,9,"choose",true)<<endl;
    cout<<"Should have 1\n";
    
    cout<<"\n========================================================\n";
    cout<<"==== SECOND TEST MYRHEIM MEYER DIMENSIONAL ESTIMATOR ====";
    cout<<"\n========================================================\n";
    Spacetime S = Spacetime();
    S.FlatSpacetime();

    cout << "\n1. CHECK MMDIM SAVING FILE FOR PYTHON\n";
    print_vector(cards);

    std::stringstream stream0;
    stream0 << std::fixed << Nreps;
    std::string Nreps_str = stream0.str();
    std::stringstream stream1;
    stream1 << std::fixed << cards[cards.size()-1];
    std::string maxsize_str = stream1.str();
    std::string filename = "../../../data/test_MMdim_forpy/MMdim_Flat_Nreps"
                            +  Nreps_str
                            + "_UpTo" + maxsize_str
                            + ".txt";

    if (run_on_hpc)
    {
      std::cout<<"Run on HPC : true\n";
      const char* homeDir = getenv("HOME");
      filename = std::string(homeDir) 
                            + "/MSci_Schwarzschild_Causets/data/test_MMdim_forpy/MMdim_Flat_Nreps"
                            +  Nreps_str
                            + "_UpTo" + maxsize_str
                            + ".txt";
    }
    std::ofstream out;
    std::cout<<"Will save in "<<std::endl;
    std::cout<<filename<<std::endl;

    // First row -> sizes
    out.open(filename);
    for (int card : cards){
      out<<card<<"avg, "<<card<<"std, ";
    }

    // Results
    for (int dim : dims)
    {
      std::vector<double> center (dim, 0.0);
      CoordinateShape shape(dim,"bicone",center,radius);
      
      out<<std::endl;
      for (int card : cards){
        std::cout<<std::endl;
        std::cout<<dim<<"D, N = "<<card<<": ";
        double sum = 0;
        double sum2 = 0;
        double nsuccess = 0;

        for (int i = 0; i<Nreps; i++){
          std::cout << i << ", ";
          SprinkledCauset Cs(card, S, shape, poisson,
                              make_matrix, special, use_transitivity,
                              make_sets, make_links, sets_type);
          int size_min = (card/5 < 100)? card/5 : 100;
          vector<double> MMd_result = Cs.MMdim_est("big", 10, 
                                                  size_min, card, true);
          
          //print(MMd_result);
          if (MMd_result[0]>0)
          {
              sum  += MMd_result[0];
              sum2 += MMd_result[0]*MMd_result[0];
              nsuccess += 1;
          }
        }

        double avg = sum/nsuccess;
        double std = std::sqrt(sum2/nsuccess - avg*avg);
        //std::cout << "Dim estimated=" <<avg <<"+-"<<std<<std::endl;
        out << avg << ", " << std;
        if (card!=cards[cards.size()-1])
        out<<",";
      }
    }
    out.close();

    auto stop = high_resolution_clock::now();
    double duration = duration_cast<microseconds>(stop - start).count();
    cout << "\n=============Time taken: "
            << duration/pow(10,6) << " seconds===============\n" << std::endl;

    std::system("python 'tests/test_causets_general/MMdim_flat_testplot.py'");

    return 0;
};