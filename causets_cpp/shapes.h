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
        std::vector<double> center;
        std::map<char, double> params;
        struct dict;
        double volume;

        CoordinateShape(int dim, char name, 
                        std::vector<double> center = {0},
                        double radius   = 1,
                        double edge     = 1,
                        std::vector<double> edges = {1},
                        double hollow   = 0,
                        double duration = 2);
        
        void param_rangecheck(char name, 
                     double maxValue = std::numeric_limits<double>::infinity(),
                     bool canBeZero = false);
        
        double              Parameter (char key);
        double              Volume();
        std::vector<double> Limits(int dimension);

};

#endif /* SHAPES_H */
