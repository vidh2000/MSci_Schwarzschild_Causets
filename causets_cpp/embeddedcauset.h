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
        vector<vector<double>> _coords;
        CoordinateShape _shape;
        Spacetime _spacetime;

        // CONSTRUCTORS
        EmbeddedCauset(const char* spacetime= "flat",
                       const char* shape    = "diamond",
                       int dim = 4,
                       vector<vector<double>> coordinates = {{0}}, 
                       bool make_matrix = true,
                       bool make_past = false,
                       bool make_fut = false,
                       bool make_past_links = false,
                       bool make_fut_links = false);
        EmbeddedCauset(const char* spacetime,
                       const char* shape,
                       int dim,
                       int size,
                       bool make_matrix = true,
                       bool make_past = false,
                       bool make_fut = false,
                       bool make_past_links = false,
                       bool make_fut_links = false);
        EmbeddedCauset(Spacetime spacetime, 
                       CoordinateShape shape, 
                       vector<vector<double>> coordinates,
                       bool make_matrix = true,
                       bool make_past = false,
                       bool make_fut = false,
                       bool make_past_links = false,
                       bool make_fut_links = false);
        EmbeddedCauset(Spacetime spacetime, 
                       CoordinateShape shape, 
                       int size,
                       bool make_matrix = true,
                       bool make_past = false,
                       bool make_fut = false,
                       bool make_past_links = false,
                       bool make_fut_links = false);
        

        // GETTERS
        int spacetime_dim();
        double density();
        double length_scale();
        CoordinateShape &Shape = _shape;
        Spacetime &Spacetime = _spacetime;
        

        // Methods of constructing the causal set
        void make_cmatrix(const char* method = "coordinates");
        void make_all_pasts(const char* method = "coordinates");
        void make_all_futures(const char* method = "coordinates");
        void make_pasts(const char* method = "coordinates");
        void make_futures(const char* method = "coordinates");
        void make_past_links(const char* method = "coordinates");
        void make_fut_links(const char* method = "coordinates");



        // RELATIONS
        bool areTimelike(vector<double> xvec, vector<double> yvec);
        bool AprecB(vector<double> xvec, vector<double> yvec);


        // MODIFY
        Causet create (vector<vector<double>> coordinates);
        void relate ();
        void relabel(const char* method = "0", bool reverse = false);      
};

#endif /*EMBEDDEDCAUSET_H*/