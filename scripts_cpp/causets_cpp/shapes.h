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

using std::vector;

class CoordinateShape
/**
 * @brief Handles a coordinate shape for embedding (and sprinkling) regions.
 */
{
    public:
        int _dim = 0;
        const char* _name;
        vector<double> _center;
        std::map<const char*, double> _params;
        double _volume = 0;

        CoordinateShape(int dim = 4,
                        const char* name = "diamond", 
                        vector<double> center = {0},
                        double radius   = 1,
                        double edge     = 1,
                        vector<double> edges = {1},
                        double hollow   = 0,
                        double duration = 2);
        
        void param_rangecheck(const char* name, 
                     double maxValue = std::numeric_limits<double>::infinity(),
                     bool canBeZero = false);
        
        double Parameter (const char* key);
        vector<double> Edges ();
        double Volume();

        vector<double> Limits(int dimension);

        // The following functions are only used in plotting, 
        // (apart for EllipseEdge in DeSitter XY slicing)
        // As these are not implemeneted, functions are left undefined:
        // double              MaxEdgeHalf();
        // vector<double> RectangleEdge();
        // vector<vector<double>> CircleEdge();
        // vector<vector<double>> ElipseEdge();
        // vector<vector<double>> BiconeEdge();
        // vector<vector<double>> CylinderCutEdge();
        // vector<vector<double>> BallSurface();
        // vector<vector<double>> CuboidSurface();
        // vector<vector<double>> CubeSurface();
        // vector<vector<double>> CylinderlSurface();
        // vector<vector<double>> OpenConeSurface();


};

#endif /* SHAPES_H */
