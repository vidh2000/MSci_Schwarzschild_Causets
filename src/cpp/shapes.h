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

/**
 * @brief Sets the embedding shape, but does not adjust event coordinates.
        The dimension parameter dim must be an integer greater than zero.

 * @param dim: int > 0.
    
 * @param name: (const char*). 
        Accepted shape names (and parameters)
        -------------------------------------
        "ball" ("radius": float, default: 1.0, must be > 0.0,
                "hollow": float, default: 0.0, must be >= 0.0 and < 1.0)
            Ball shape in all spacetime coordinates. 
        "bicone" or"diamond"("radius": float, default: 1.0, must be > 0.0,
                  "hollow": float, default: 0.0, must be >= 0.0 and < 1.0)
            Ball shape in all space coordinates and conical to the past and 
            future. 
        "cylinder" ("radius": float, default: 1.0, must be > 0.0,
                    "duration": float, default: 2.0 * radius, must be > 0.0,
                    "hollow": float, default: 0.0, must be >= 0.0 and < 1.0)
            Ball shape in all space coordinates and straight along the time 
            coordinate for the length "duration". 
        "cube" ("edge": float, default: 1.0, must be > 0.0)
            Cube shape with the same edge length "edge" in all spacetime 
            coordinates.
        "cuboid" ("edges": Iterable[float], default: [1.0, 1.0, ...], 
                                            must all be > 0.0)
            Cuboid shape with distinct edge lengths "edges" in the respective 
            spacetime coordinates. The default edges yield a cube.
    
    @param center: vector of coordinates of center. Default 0\vec.

    @param radius: float > 0 for "ball", "bicone", "cylinder". Default 1.

    @param duration: time extension of cylinder. Default 2.

    @param hollow: fraction [0,1) of interior which is hollow. Default 0.

    @param edge: float > 0 for "cube". Default 1.

    @param edges: vector of floats > 0 for "cuboid". Default {1}.

 * 
 */
class CoordinateShape
/**
 * @brief Handles a coordinate shape for embedding (and sprinkling) regions.
 */
{
    public:
        int _dim;
        const char* _name;
        std::vector<double> _center;
        std::map<std::string, double> _params;
        double _volume = 0;
        double _deltaphi= 0;
        bool isBicone;

        //CONSTRUCTOR (& RELATED)
        CoordinateShape(int dim = 4,
                        const char* name = "bicone", 
                        std::vector<double> center = {},//{0}
                        double radius   = 1,
                        double duration = 2,
                        double hollow   = 0,
                        double edge     = 1,
                        std::vector<double> edges = {});//{1}
        

        void limitphi(double deltaphi);
                        
                        
        //INTERNAL WORKINGS

        void param_rangecheck(std::string name, 
                              double maxValue = std::nan(""),
                              bool canBeZero = false);

        
        
        //GETTERS

        double Parameter (std::string key);
        std::vector<double> Edges ();
        double Volume();
        std::vector<double> Limits(int dimension);
};

#endif /* SHAPES_H */
