/// \authors Vid Homsak, Stefano Veroni
/// \date 31/10/2022

#ifndef CAUSET_H
#define CAUSET_H

#include <cmath>
#include <set>
#include <vector>
#include <unordered_set>
#include <stdint.h>


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
        std::vector<std::vector<int>> _CMatrix = {};
        bool _special_matrix = false;
        int _size = 0;
        int _dim = 0;
        std::vector<std::unordered_set<int>> _pasts   = {};
        std::vector<std::unordered_set<int>> _futures = {};
        std::vector<std::unordered_set<int>> _past_links   = {};
        std::vector<std::unordered_set<int>> _future_links = {};


        // CONSTRUCTOR

        Causet(std::vector<std::vector<int>> Cmatrix = {});
        
        template <typename num>
        Causet(std::vector<std::vector<num>> Cmatrix);
        

        //SETTERS/GETTERS

        std::vector<std::vector<int>> get_CMatrix();
        
        // MAKERS

        void make_sets_fromC();
        virtual void make_cmatrix();
        virtual void make_pasts();
        virtual void make_futures();
        virtual void make_past_links();
        virtual void make_future_links_fromC();

        std::vector<std::vector<int>> CMatrix(std::vector<int> labels = {});
        int size();
        bool is_CMatrix_special();
        bool is_Cij_special();


        // RELATIONS (not supported anymore I think as trivial)

        bool areTimelike(int a, int b);
        bool AprecB(int a, int b);


        // KINEMATICS

        static double ord_fr(Causet & A,
                      const char* denominator = "choose");
        static double ord_fr(std::vector<std::vector<int>> & A,
                      const char* denominator = "choose");
        template<typename SET>
        static double ord_fr(std::vector<SET> & A_pasts,
                      const char* denominator = "choose");
        double ord_fr(int a, int b,
                      const char* denominator = "choose",
                      bool from_matrix = true);

        static double MM_drelation(double d); 
        
        std::vector<double> MMdim_est(const char* method = "random",
                        int Nsamples = 20,
                        int size_min = 100,
                        double size_max = 1e9,
                        bool from_matrix = true);   

        std::vector<std::vector<int>> getIntervalCmatrix(
                                std::vector<int> ordered_interval);
                                
        std::vector<double> Nk_BD (int x, int kmax, int kmin = 1);
        


        // INTERVAL

        Causet Interval(int a, int b,
                        bool includeBoundary = true,
                        bool disjoin = false,
                        const char* createmethod = "set");
        int IntervalCard(int a, int b, bool includeBoundary = true);


        // CAUSET REPRESENTATION & SAVING (to be added...)

        std::vector<std::vector<double>> LMatrix (std::vector<int> labels = {});
        void saveC(const char* path_file_ext);


        // MODIFIERS

        void coarsegrain(int card, bool make_matrix = true, 
                     bool make_sets = false, bool make_links = true,
                     int seed = 0);
        void cgrain(int card, bool make_matrix = true, 
                     bool make_sets = false, bool make_links = true,
                     int seed = 0);
        void coarsegrain(double fract, bool make_matrix = true, 
                     bool make_sets = false, bool make_links = true,
                     int seed = 0);
        void cgrain(double fract, bool make_matrix = true, 
                     bool make_sets = false, bool make_links = true,
                     int seed = 0);

        void discard(int label, bool make_matrix = true, 
                     bool make_sets = false, bool make_links = true);  
        void discard(std::vector<int> labels, bool make_matrix = true, 
                     bool make_sets = false, bool make_links = true);  

        //Destructor
        //virtual ~Causet();

        //MAKE ATTRS


};

#endif /*CAUSET_H*/