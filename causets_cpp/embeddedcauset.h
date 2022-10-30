/// \authors Vid Homsak, Stefano Veroni
/// \date 29/09/2022

#ifndef EMBEDDEDCAUSET_H
#define EMBEDDEDCAUSET_H

#include <set>
#include <vector>

#include "causet.h"
#include "spacetimes.h"
#include "shapes.h"

using std::vector;
using std::set;

/**
 * @brief An embedded causet in a spacetime subset of a specified shape.
 */
class EmbeddedCauset: public Causet
{
    public:
        CoordinateShape _shape;
        Spacetime _spacetime;

        // CONSTRUCTORS
        EmbeddedCauset(const char* spacetime= "flat",
                       const char* shape    = "diamond",
                       int dim = 4,
                       vector<vector<double>> coordinates = {{0}});
        EmbeddedCauset(const char* spacetime,
                       CoordinateShape shape,
                       int dim,
                       vector<vector<double>> coordinates);
        EmbeddedCauset(Spacetime spacetime,
                       const char* shape, 
                       vector<vector<double>> coordinates);
        EmbeddedCauset(Spacetime spacetime, 
                       CoordinateShape shape, 
                       vector<vector<double>> coordinates);
        
        // GETTERS
        int Dim();
        double Density();
        double LengthScale();
        CoordinateShape &Shape = _shape;
        Spacetime &Spacetime = _spacetime;
        
        // MODIFY
        Causet create (vector<vector<double>> coordinates);
        void relate ();
        void relabel(int guiding_dim = 0, bool reverse = false);      
};

#endif /*EMBEDDEDCAUSET_H*/