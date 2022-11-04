/// \authors Vid Homsak, Stefano Veroni
/// \date 28/09/2022

#ifndef SPRINKLEDCAUSET_H
#define SPRINKLEDCAUSET_H

#include <set>
#include <vector>

#include "causet.h"
#include "embeddedcauset.h"
#include "spacetimes.h"
#include "shapes.h"

using std::vector;
using std::set;

/**
 * @brief An causet obtained from sprinkling in a spacetime subset 
 *        of a specified shape.
 */
class SprinkledCauset: public EmbeddedCauset
{
    public:
        //Inherited
        // // vector<vector<int8_t>> _CMatrix;
        // // int _size = 0;
        // // int _dim = 0;
        // // bool _special_matrix;
        
        // // vector<std::unordered_set<int>> _pasts   = {};
        // // vector<std::unordered_set<int>> _futures = {};
        // // vector<std::unordered_set<int>> _past_links   = {};
        // // vector<std::unordered_set<int>> _future_links = {};
        // vector<vector<double>> _coords;
        // CoordinateShape _shape;
        // Spacetime _spacetime;
        double _intensity = -1;
       
        // Constructor
        SprinkledCauset(int card,
                        Spacetime spacetime, 
                        CoordinateShape shape, 
                        bool poisson = false,
                        bool make_matrix = true,
                        bool special = false,
                        bool use_transitivity = true,
                        bool make_sets = false,
                        bool make_links = false,
                        const char* sets_type = "past",
                        int seed = std::nanf(""));
        
        // Properties
        double &Intensity = _intensity;

        // Methods
        static vector<vector<double>> sprinkle_coords(int count,
                                                 CoordinateShape shape,
                                                 int seed = std::nanf(""));
        static vector<vector<double>> sprinkle(int count, CoordinateShape shape, 
                                                bool poisson = false,
                                                int seed = std::nanf(""));
     
};

#endif /* SPRINKLEDCAUSET_H */