#include <algorithm>
#include <cmath>
#include <cstdio>
#include <chrono>
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

#include "../causets_cpp/sprinkledcauset.h"
#include "../causets_cpp/shapes.h"
#include "../causets_cpp/spacetimes.h"

#include "../causets_cpp/functions.h"
#include "../causets_cpp/vecfunctions.h"

#include "../causets_cpp/kinematics_functions.h"

//#include <boost/range/combine.hpp>
//#include <boost/math/special_functions/gamma.hpp>
#include <omp.h>

// $HOME var get
#include <unistd.h>
#include <sys/types.h>
//#include <pwd.h>

using namespace std::chrono;
using namespace boost::math;





// Defining Ricci scalar (components) and proper time extension of the interval
///////////////////////////////////////////////////////////////////////////////





int main(){


// Create causet in Schwarzschild spacetime
///////////////////////////////////////////////////////////////////////////////
                
std::vector<int> dims = {4}; 
std::vector<int> cards = {1000};
int min_size = 5;  //Minimal size of the interval (min # of elements in it)
int max_size = 0;   //if 0 -> ==_size of the causet
double mass = 0.25;
int N_reps = 5;
int N_intervals = 10; // Number of intervals per causet realisation


// Sprinkling Parameters
///////////////////////////////////////////
bool poisson          = true;
bool make_matrix      = true;
bool special          = false;
bool use_transitivity = false;
bool make_sets        = false;
bool make_links       = false; 
const char* sets_type = "future";
const char* name      = "cylinder";

// Shape parameters
double radius = 5*mass;
double height = 1; 

///////////////////////////////////////////

// Begin program
auto beginning = high_resolution_clock::now();

std::cout<<"\n\n============= Sprinkling into "<<name<<" ==================\n";
int iteration = 0;
for (auto dim: dims)
{
    auto start = high_resolution_clock::now();

    for (auto card : cards)
    {
        /*Vectors where to store results*/
        std::vector<double> vec_of_RicciScalars (N_reps*N_intervals);
        std::vector<double> vec_of_Ricci00s     (N_reps*N_intervals); 
        std::vector<double> vec_of_radii        (N_reps*N_intervals);
        std::vector<double> vec_of_RSSMMdim     (N_reps);
        std::vector<double> vec_of_RSSMMstd     (N_reps);

        std::vector<double> vec_of_RicciBD      (N_reps*N_intervals);
        std::vector<double> vec_of_radiiBD      (N_reps*N_intervals);

         // Set up shape for repetitions
        std::vector<double> center (dim, 0.);
        CoordinateShape shape(dim,name,center,radius,height);

        // Set up spacetime for repetitions
        Spacetime S = Spacetime();
        S.BlackHoleSpacetime(dim,mass);

        // Generate causets and find values
        ////////////////////////////////////////////////////////////////////
        int rep=0;
        while (rep<N_reps)
        {
            std::cout << "Dim="<< dim 
                      <<", Card="<< card<<", "
                      <<"Minsize="<<min_size<<", "
                      <<"Rep: "<<(rep+1)<<"/"<<N_reps<<"\n";
            auto repstart = high_resolution_clock::now();

            // Sprinkle the causet
            SprinkledCauset C(card, S, shape, poisson,
                              make_matrix, special,    use_transitivity,
                              make_sets,   make_links, sets_type);
            
            // Get result - N_intervals RicciScalars, Ricci00s, radii - 
            // of nth repetition, with RSS model
            std::vector<std::vector<double>> nresults;
            std::vector<double> RSSMMdim;
            try {
                std::cout<<"\nDoing RSS now";
                nresults = Riccis_RSS_radial_sample(C, N_intervals,
                                                     min_size, max_size);
                RSSMMdim = RSSMMdim_sample(C, N_intervals,
                                                     min_size, max_size);
            } catch (std::runtime_error& e) {
                continue;
            }  

            // and with BD model
            std::cout<<"\nDoing BD now";
            std::vector<std::vector<double>> nBDresults 
            = R_BD_sample(C, N_intervals);

            // Add to all results
            std::cout<<"\nChecking parallelisation and random works: radii"<<
                        " from BD are : ";
            for (int i = 0; i<N_intervals; i++)
            {
                vec_of_RicciScalars[rep*N_intervals+i] = nresults[i][0];
                vec_of_Ricci00s    [rep*N_intervals+i] = nresults[i][1];
                vec_of_radii       [rep*N_intervals+i] = nresults[i][2];

                vec_of_RicciBD     [rep*N_intervals+i] = nBDresults[i][0];
                vec_of_radiiBD     [rep*N_intervals+i] = nBDresults[i][1];
                std::cout<<nBDresults[i][1]<< ", ";
            }
            vec_of_RSSMMdim[rep] = RSSMMdim[0];
            vec_of_RSSMMstd[rep] = RSSMMdim[1];

            //Timing rep
            auto repend = high_resolution_clock::now();
            double duration = duration_cast<microseconds>(repend - repstart).count();
            std::cout << "Rep finished, N = " << C._size
            << ": " << duration/pow(10,6) << " seconds" << std::endl;

            // Add to counter of repetitions
            rep++;
        } /*end of loop over different relaisations*/
    
        auto mid = high_resolution_clock::now();
        double duration = duration_cast<microseconds>(mid - start).count();
        
        std::cout << "-------------------------------------\n";
        std::cout << "<t> for "<< N_reps
                << " causets in D=" << dim << " with N = " << card << ": "  
                << duration/pow(10,6)/N_reps
                << " seconds\n" << std::endl;

        //////////////////////////////////////////////////////////////////////
        ///// ANALYSIS OF RESULTS
        //////////////////////////////////////////////////////////////////////

        // Calculate the density - assume we always sprinkle into cylinder
        double density, volume;
        if (dim==3) volume = pi*radius*radius*height;
        else if (dim==4) volume = 4/3*pi*radius*radius*radius*height;
        else {
            print("ERROR: Dimension must be 3 or 4!");
            throw std::runtime_error("");
        }
        density = card/volume;

        
        print("\nxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n");
        std::cout << "Density = " << density << std::endl;
        print("Using units of discreteness scale");

        // get number of bins for Friedman-Diaconis rule
        std::vector<double> sorted_rs = vec_of_radii;
        sort(sorted_rs.begin(), sorted_rs.end());
        int N = sorted_rs.size();
        double q1 = sorted_rs[N/4];
        double q3 = sorted_rs[3*N/4];
        double iqr = q3 - q1;
        double bin_width = 2 * iqr * pow((double)N, -1.0/3.0);
        int n_rad_bins = (sorted_rs[N] - sorted_rs[0])/bin_width;
        std::cout<<"FD rule"<<std::endl;

        // get average of RSSMMdim
        std::vector<int> Ns (N_reps, (int)N_intervals);
        std::pair<double, double>  RSSMMdimfinal = combine_meass(
                                Ns, vec_of_RSSMMdim, vec_of_RSSMMstd);
        double MMdim_avg = RSSMMdimfinal.first;
        double MMdim_est = RSSMMdimfinal.second;
        std::cout<<"\nDIMENSION ESTIMATES IS :"<<MMdim_avg<<" +- "<<MMdim_est
                 <<std::endl;

        // get average R and R00 per bin
        std::vector<double> avg_R (n_rad_bins);
        std::vector<double> std_R (n_rad_bins);
        std::vector<double> avg_R00 (n_rad_bins);
        std::vector<double> std_R00 (n_rad_bins);

        std::vector<double> avg_RBD (n_rad_bins);
        std::vector<double> std_RBD (n_rad_bins);

        std::vector<int>    bin_counts (n_rad_bins);
        std::cout<<"\nbins first start";
        for (int i = 0; i < N; i++)
        {
            int bin = (vec_of_radii[i] - sorted_rs[0])/bin_width;
            std::cout<<bin<<std::endl;
            avg_R[bin]      += vec_of_RicciScalars[i];
            std_R[bin]      += vec_of_RicciScalars[i]*vec_of_RicciScalars[i];
            avg_R00[bin]    += vec_of_Ricci00s[i];
            std_R00[bin]    += vec_of_Ricci00s[i]*vec_of_Ricci00s[i];

            avg_RBD[bin]      += vec_of_RicciBD[i];
            std_RBD[bin]      += vec_of_RicciBD[i]*vec_of_RicciBD[i];

            bin_counts[bin] += 1;
        }
        std::cout<<"\nbins first";
        printf( "| %11s | %19s | %19s | %19s \n", 
            "Bin","Ricci Scalar (RSS)","Ricci Scalar (BD)","Ricci 00 (RSS)");
        for (int bin = 0; bin < n_rad_bins; bin++)
        {
            float r1 = sorted_rs[0] + bin*bin_width;
            float r2 = sorted_rs[0] + (bin+1)*bin_width;
            avg_R[bin]      /= bin_counts[bin];
            std_R[bin]      /= bin_counts[bin];
            std_R[bin]      -= avg_R[bin]*avg_R[bin];
            std_R[bin]       = std::pow(std_R[bin], 0.5);
            avg_R00[bin]    /= bin_counts[bin];
            std_R00[bin]    /= bin_counts[bin];
            std_R00[bin]    -= avg_R00[bin]*avg_R00[bin];
            std_R00[bin]     = std::pow(std_R00[bin], 0.5);

            avg_RBD[bin]    /= bin_counts[bin];
            std_RBD[bin]    /= bin_counts[bin];
            std_RBD[bin]    -= avg_RBD[bin]*avg_RBD[bin];
            std_RBD[bin]     = std::pow(std_RBD[bin], 0.5);

            printf("| %4f - %4f |  %7f +- %6f  |  %7f +- %6f  |  %7f +- %6f  |\n", 
            r1, r2, avg_R[bin]  , std_R[bin], 
                    avg_RBD[bin], std_RBD[bin], 
                    avg_R00[bin], std_R00[bin]);
        }  
    } 
}

auto finish = high_resolution_clock::now();
double duration = duration_cast<microseconds>(finish - beginning).count();
std::cout<<"===============================================================\n";
std::cout << "\nProgram took in total: "
        << duration/pow(10,6) << " seconds\n" << std::endl;




// end of main   
}