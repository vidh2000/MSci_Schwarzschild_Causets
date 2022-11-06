/// \authors Vid Homsak, Stefano Veroni
/// \date 31/10/2022

#ifndef CAUSET_H
#define CAUSET_H

#include <cmath>
#include <set>
#include <vector>
#include <unordered_set>

using std::vector;
using std::set;
using std::unordered_set;

/**
 * @brief Causet class.
 *
 * 
 * @param CMatrix: a vector of vectors of integers.
 *  Essentially it is the causal matrix where C_{i,j}. The matrix can be
 * "standard" or "special" (see special_matrix attribute).
 *  If standard: C_ij = 1 if e_j<e_i, 0 otherwise.
 *  If special: C_ij = -1 if e_j<*e_i (link), 1 if e_j<<e_i(not_link), 
 *  0 otherwise.
 * @param special_matrix: true if CMatrix is special, 0 if standard.
 * @param size: the number of elements
 * @param dim: the Myrheim Meyer dimension.
 */


class Causet
{
    public:
        // Attributes
        vector<vector<int8_t>> _CMatrix = {};
        bool _special_matrix = false;
        int _size = 0;
        int _dim = 0;
        
        vector<std::unordered_set<int>> _pasts   = {};
        vector<std::unordered_set<int>> _futures = {};
        vector<std::unordered_set<int>> _past_links   = {};
        vector<std::unordered_set<int>> _future_links = {};


        // CONSTRUCTOR
        Causet ();
        Causet(vector<vector<double>> Cmatrix, 
                bool past_links = false, bool fut_links = false);
        Causet(vector<vector<double>> coordinates,
               const char* method = "pasts");
        

        //SETTERS/GETTERS
        void make_cmatrix();
        void make_lmatrix();
        void make_pasts();
        void make_futures();
        void make_past_links();
        void make_future_links();

        vector<vector<int>> CMatrix(set<int> labels = {});
        int size();
        bool is_CMatrix_special();
        bool is_Cij_special();


        // RELATIONS
        bool areTimelike(int a, int b);
        bool AprecB(int a, int b);


        // KINEMATICS
        static double ord_fr(Causet A,
                      const char* denominator = "choose",
                      bool isdisjoined = true);
        static double ord_fr(vector<vector<int8_t>> A,
                      const char* denominator = "choose",
                      bool isdisjoined = true);
        template<typename SET>
        static double ord_fr(vector<SET> A_futures, 
                      vector<SET> A_pasts,
                      const char* denominator = "choose",
                      bool isdisjoined = true);
        double ord_fr(int a, int b,
                      const char* denominator = "choose");

        static double MM_drelation(double d); 
        
        vector<double> MMdim_est(const char* method = "random",
                        int d0 = 2, int Nsamples = 20,
                        int size_min = 10, double size_max = nan(""));   
        


        // INTERVAL
        Causet Interval(int a, int b,
                        bool includeBoundary = true,
                        bool disjoin = false,
                        const char* createmethod = "set");
        int IntervalCard(int a, int b, bool includeBoundary = true);


        // CAUSET REPRESENTATION & SAVING (to be added...)
        vector<vector<double>> LMatrix (set<int> labels = {});
        void saveCasCSV(const char* filename);
        void saveCasTXT(const char* filename);
        void saveLasCSV(const char* filename);
        void saveLasTXT(const char* filename);


        // MODIFIERS
        void coarsegrain(int card, bool make_matrix = true, 
                     bool make_sets = false, bool make_links = true);
        void cgrain(int card, bool make_matrix = true, 
                     bool make_sets = false, bool make_links = true);
        void coarsegrain(double fract, bool make_matrix = true, 
                     bool make_sets = false, bool make_links = true);
        void cgrain(double fract, bool make_matrix = true, 
                     bool make_sets = false, bool make_links = true);

        void discard(int label, bool make_matrix = true, 
                     bool make_sets = false, bool make_links = true);  
        void discard(vector<int> labels, bool make_matrix = true, 
                     bool make_sets = false, bool make_links = true);
        
        static void discard_from_set(unordered_set<int> &myset, int label);
        template<typename m>
        static void discard_from_set(unordered_set<m> &myset, 
                                        vector<int> labels);
        
        static vector<int> distinct_int_random (int card, int N, int seed = 0);
        static vector<int> distinct_int_random1(int card, int N, int seed = 0);
        static vector<int> distinct_int_random2(int card, int N, int seed = 0);
        


};

#endif /*CAUSET_H*/