#ifndef SHAPES_H
#define SHAPES_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>
#include <fstream>
#include <stack>
#include <string>
#include <stdio.h>
#include <vector>


class CoordinateShape
/**
 * @brief Handles a coordinate shape for embedding (and sprinkling) regions.
 */
{
    public:
        int dim;
        char name;
        std::vector<float> center;
        std::map<char, float> params;
        struct dict;
        float volume;

        CoordinateShape(int dim, char name, 
                        std::vector<float> center = {0},
                        float radius   = 1,
                        float hollow   = 0,
                        float duration = 2,
                        float edge     = 1,
                        std::vector<float> edges = {1});
        
        void param_rangecheck(char name, 
                     float maxValue = std::numeric_limits<double>::infinity(),
                     bool canBeZero = false);
        
        float              Parameter (char key);
        float              Volume();
        std::vector<float> Limits(int dimension);

};

#endif /* SHAPES_H */
