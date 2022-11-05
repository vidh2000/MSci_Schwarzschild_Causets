/// \authors Vid Homsak, Stefano Veroni
/// \date 29/10/2022

#ifndef EMBEDDEDCAUSET_H
#define EMBEDDEDCAUSET_H

#include <cmath>
#include <set>
#include <vector>
#include <unordered_set>

#include "causet.h"
#include "spacetimes.h"
#include "shapes.h"

using std::vector;
using std::set;
using std::unordered_set;

/**
 * @brief An embedded causet in a spacetime subset of a specified shape.
 */
class EmbeddedCauset: public Causet
{
    public:
        //Inherited
        // vector<vector<int8_t>> _CMatrix = {};
        // bool _special_matrix = false;
        // int _size = 0;
        // int _dim = 0;
        
        // vector<std::unordered_set<int>> _pasts   = {};
        // vector<std::unordered_set<int>> _futures = {};
        // vector<std::unordered_set<int>> _past_links   = {};
        // vector<std::unordered_set<int>> _future_links = {};
        vector<vector<double>> _coords;
        CoordinateShape _shape;
        Spacetime _spacetime;

        // CONSTRUCTORS
        EmbeddedCauset();
        EmbeddedCauset(Spacetime spacetime, 
                       CoordinateShape shape, 
                       vector<vector<double>> coordinates,
                       bool make_matrix = true,
                       bool special = false,
                       bool use_transitivity = true,
                       bool make_sets = false,
                       bool make_links = false,
                       const char* sets_type = "past");
        

        // Methods of constructing the causal set attributes
        void EmbeddedCauset::make_cmatrix (const char* method = "coordinates",
                                    bool special = false,
                                    bool use_transitivity = true,
                                    bool make_sets = false,
                                    bool make_links = false,
                                    const char* sets_type = "past");
        void make_all_pasts  (const char* method = "coordinates");
        void make_all_futures(const char* method = "coordinates");
        void make_pasts      (const char* method = "coordinates");
        void make_futures    (const char* method = "coordinates");
        void make_past_links (const char* method = "coordinates");
        void make_fut_links  (const char* method = "coordinates");


        // GETTERS
        int spacetime_dim();
        double density();
        double length_scale();
        CoordinateShape &get_shape = _shape;
        Spacetime &get_spacetime = _spacetime;


        // RELATIONS
        bool areTimelike(vector<double> xvec, vector<double> yvec);
        bool AprecB(vector<double> xvec, vector<double> yvec);


        // MODIFIERS
        void sort_coords(int dim = 0, bool reverse = false);
        void relabel(const char* method = "0", bool reverse = false);   
        void add(vector<double> xvec);
        void discard(int label, bool make_matrix = true, 
                     bool make_sets = false, bool make_links = true);  
        void discard(vector<int> labels, bool make_matrix = true, 
                     bool make_sets = false, bool make_links = true);
        
        template<typename n>
        static void discard_from_set(unordered_set<n> &myset,
                                     n label);
        template<typename n>
        static void discard_from_set(unordered_set<n> &myset, 
                                     vector<int> labels);
};

#endif /*EMBEDDEDCAUSET_H*/