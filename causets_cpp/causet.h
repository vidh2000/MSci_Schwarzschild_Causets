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
        vector<vector<int8_t>> _CMatrix;
        int _size = 0;
        int _dim = 0;
        bool _special_matrix;
        
        // Creating useful "representations" of causet
        //vector<std::unordered_set<int>> pasts;
        vector<std::unordered_set<int>> _pasts   = {};
        vector<std::unordered_set<int>> _futures = {};
        vector<std::unordered_set<int>> _past_links   = {};
        vector<std::unordered_set<int>> _future_links = {};

        // CONSTRUCTOR
        Causet(vector<vector<double>> Cmatrix, 
                bool past_links = false, bool fut_links = false);
        Causet(vector<vector<double>> coordinates,
               const char* method = "pasts");
        
        //SETTERS/GETTERS
        void make_Cmatrix();
        void make_Lmatrix();
        void make_pasts();
        void make_futures();
        void make_past_links();
        void make_future_links();

        vector<vector<int>> CMatrix(set<int> labels = {});
        int size();
        int dim();
        bool is_CMatrix_special();
        bool is_Cij_special();


        // RELATIONS
        bool areTimelike(int a, int b);
        bool AprecB(int a, int b);


        // KINEMATICS
        double ord_fr(Causet A,
                      const char* denominator = "choose",
                      bool isdisjoined = true);
        double ord_fr(vector<vector<int8_t>> A,
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
        // typedef double (*func)();//need to pick optimiser;
        template <typename F>
        double MMdim_est(const char* method = "random",
                        int d0 = 2, int Nsamples = 20,
                        int size_min = 10, double size_max = nan(""),
                        F optimiser = Causet::optimiser_placeholder);   

        // INTERVAL
        Causet Interval(int a, int b,
                        bool includeBoundary = true,
                        bool disjoin = false,
                        const char* createmethod = "set");
        int IntervalCard(int a, int b, bool includeBoundary = true);

        // CAUSET REPRESENTATION & SAVING (to be added...)
        void relabel();
        vector<vector<double>> LMatrix (set<int> labels = {});
        void saveCasCSV(const char* filename);
        void saveCasTXT(const char* filename);
        void saveLasCSV(const char* filename);
        void saveLasTXT(const char* filename);

};

#endif /*CAUSET_H*/