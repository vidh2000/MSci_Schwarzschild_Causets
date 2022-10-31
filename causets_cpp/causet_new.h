/// \authors Vid Homsak, Stefano Veroni
/// \date 31/10/2022

#ifndef CAUSET_H
#define CAUSET_H

#include <cmath>
#include <set>
#include <vector>

using std::vector;
using std::set;

/**
 * @brief Causet class.
 *
 * 
 * @param causet: a vector of vectors of integers.
 *  Essentially it is the causal matrix where (M)_{i,j}
 *  can take value of 1 if e_j<e_i, 0 if they aren't related and
 *  -1 if e_j<e_i and they are a link. //not implemented
 */

class Causet
{
    public:
        
        // Attributes
        vector<vector<double>> coords;
        vector<vector<int>> causet;
        int size;
        int dim;
        
        // Creating useful "representations" of causet
        vector<set<int>> pasts;
        vector<set<int>> futures;
        vector<set<int>> past_links;
        vector<set<int>> future_links;

        // CONSTRUCTOR
        Causet(vector<vector<double>> coordinates,
               const char* method = "pasts");
        
        // Methods of constructing the causal set
        void make_pasts();
        void make_cmatrix();

        // KINEMATICS
        double ord_fr(Causet A,
                      const char* denominator = "choose",
                      bool isdisjoined = true);
        double ord_fr(vector<set<int>> A_future = {}, 
                      vector<set<int>> A_past = {},
                      const char* denominator = "choose",
                      bool isdisjoined = true);
        double ord_fr(int a, int b,
                      const char* denominator = "choose",
                      bool isdisjoined = true);
        static double optimiser_placeholder();
        typedef double (*func)();//need to pick optimiser;
        double MMdim_est(const char* method = "random",
                        int d0 = 2, int Nsamples = 20,
                        int size_min = 10, double size_max = nan(""),
                        func optimiser = Causet::optimiser_placeholder);   

        // INTERVAL
        Causet Interval(int a, int b,
                        bool includeBoundary = true,
                        bool disjoin = false,
                        const char* createmethod = "set");
        int IntervalCard(int a, int b, bool includeBoundary = true);

        // CAUSET REPRESENTATION & SAVING (to be added...)
        vector<vector<double>> CMatrix (const char* method="causality",
                                         set<int> labels = {});
        vector<vector<double>> CTMatrix (const char* method="causality",
                                         set<int> labels = {});
        vector<vector<double>> LMatrix (const char* method="causality",
                                         set<int> labels = {});
        void saveCasCSV(const char* filename);
        void saveCasTXT(const char* filename);
        void saveLasCSV(const char* filename);
        void saveLasTXT(const char* filename);

};

#endif /*CAUSET_H*/