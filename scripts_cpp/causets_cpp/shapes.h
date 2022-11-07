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
        int _dim = 0;
        const char* _name;
        std::vector<double> _center;
        std::map<const char*, double> _params;
        double _volume = 0;

        CoordinateShape(int dim = 4,
                        const char* name = "diamond", 
                        std::vector<double> center = {0},
                        double radius   = 1,
                        double edge     = 1,
                        std::vector<double> edges = {1},
                        double hollow   = 0,
                        double duration = 2);
        
        void param_rangecheck(const char* name, 
                                double maxValue = std::nan(""),
                                bool canBeZero = false);
        
        double Parameter (const char* key);
        std::vector<double> Edges ();
        double Volume();

        std::vector<double> Limits(int dimension);

        // The following functions are only used in plotting, 
        // (apart for EllipseEdge in DeSitter XY slicing)
        // As these are not implemeneted, functions are left undefined:
        // double              MaxEdgeHalf();
        // std::vector<double> RectangleEdge();
        // std::vector<std::vector<double>> CircleEdge();
        // std::vector<std::vector<double>> ElipseEdge();
        // std::vector<std::vector<double>> BiconeEdge();
        // std::vector<std::vector<double>> CylinderCutEdge();
        // std::vector<std::vector<double>> BallSurface();
        // std::vector<std::vector<double>> CuboidSurface();
        // std::vector<std::vector<double>> CubeSurface();
        // std::vector<std::vector<double>> CylinderlSurface();
        // std::vector<std::vector<double>> OpenConeSurface();


};

#endif /* SHAPES_H */
