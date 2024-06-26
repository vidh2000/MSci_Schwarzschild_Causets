/// \authors Vid Homsak, Stefano Veroni
/// \date 28/09/2022

#ifndef SPRINKLEDCAUSET_H
#define SPRINKLEDCAUSET_H

#include <vector>
#include "causet.h"
#include "embeddedcauset.h"
#include "spacetimes.h"
#include "shapes.h"

/**
 * @brief Construct a new Sprinkled Causet object. 
 * 
 * @param card: int, number of events to sprinkle
 * @param spacetime: Spacetime.
 * @param shape: CoordinateShape.
 * @param poisson: bool, if true card is used as average of a Poisson
 * distribution to get a new cardinality, otherwise card is used.
 * @param make_matrix: bool, if true make matrix.
 * @param special: bool, if true (and make_matrix) make special Cmatrix
 * @param use_transitivity: bool, if true use also transitivity to establish
 * causality relations. 
 * @param make_sets: bool, if true make set (see sets_type)
 * @param make_links: bool, if true make links (see sets_type)
 * @param sets_type: const char* specifying the type of set:
 * - "past": make _past_links
 * - "future": make _future_links
 * @param seed: int, for coordinate generation.
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
       
        // CONSTRUCTOR
        SprinkledCauset(double card,
                        Spacetime spacetime, 
                        CoordinateShape shape, 
                        bool poisson = true,
                        bool make_matrix = true,
                        bool special = false,
                        bool use_transitivity = true,
                        bool make_sets = false,
                        bool make_links = true,
                        const char* sets_type = "future",
                        int seed = 0,
                        int poisson_seed = 0);

        SprinkledCauset(int card,
                        Spacetime spacetime, 
                        CoordinateShape shape, 
                        bool poisson = true,
                        bool make_matrix = true,
                        bool special = false,
                        bool use_transitivity = true,
                        bool make_sets = false,
                        bool make_links = true,
                        const char* sets_type = "future",
                        int seed = 0,
                        int poisson_seed = 0);


        // METHODS
        static std::vector<std::vector<double>> sprinkle(int card, 
                                                        CoordinateShape shape, 
                                                        bool poisson = true,
                                                        int seed = 0,
                                                        int poisson_seed = 0);
        static std::vector<std::vector<double>> sprinkle(double card, 
                                                        CoordinateShape shape, 
                                                        bool poisson = true,
                                                        int seed = 0,
                                                        int poisson_seed = 0);
        static std::vector<std::vector<double>> sprinkle_coords(int count,
                                                        CoordinateShape shape,
                                                        int seed = 0);
};

#endif /* SPRINKLEDCAUSET_H */