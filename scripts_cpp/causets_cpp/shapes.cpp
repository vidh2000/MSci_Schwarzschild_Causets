#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <limits>
#include <map>
#include <math.h>
#include <numeric>
#include <stack>
#include <string>
#include <stdio.h>
#include <bits/stdc++.h>
#include <stdexcept>
#include <vector>

#include "vecfunctions.h"
//#include "functions.h"

#include "shapes.h"

#define M_PI           3.14159265358979323846  /* pi */

using std::vector;


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

    @param edge: float > 0 for "cube". Default 1.

    @param edges: vector of floats > 0 for "cuboid". Default {1}.

    @param hollow: fraction [0,1) of interior which is hollow. Default 0.

    @param duration: time extension of cylinder. Default 2.
 * 
 */
CoordinateShape::CoordinateShape(int dim, const char* name, 
                                vector<double> center,
                                double radius,
                                double edge,
                                vector<double> edges,
                                double hollow,
                                double duration)
{
    // Set Dimension
    if (dim < 1)
       {throw std::invalid_argument("Dim smaller than 1!");}
    _dim = dim;

    // Set Name
    name = (name == "diamond") ? "bicone" : name;
    if (name!="ball" || name!="bicone" || name!="cylinder"
     || name!="cube" || name!="cuboid")
     { 
        std::string dim_error = "The given shape is ";
        dim_error += name; 
        throw std::invalid_argument(dim_error + " and not supported!");
     }
    _name = name;

    // Set Center
    vector<double> zero_center (_dim, 0.0);
    if (center.size() == dim)
        {_center = center;}
    else if (center.size() == 0)
        {_center = zero_center;}
    else 
        {throw std::invalid_argument("Center's size neither 0 nor _dim");}
    
    
    // Set Shape Parameters
    if (name == "ball" or name == "bicone")
    {
        _params.insert({"radius", radius}); 
        this->param_rangecheck("radius");

        _params.insert({"hollow", hollow}); 
        this->param_rangecheck("hollow", 1.0, true);
    }
    else if (name == "cylinder")
    {
        _params.insert({"radius", radius}); 
        this->param_rangecheck("radius");

        _params.insert({"hollow", hollow}); 
        this->param_rangecheck("hollow", 1.0, true); 

        _params.insert({"duration", duration}); 
        this->param_rangecheck("duration");   
    }
    else if (name == "cube")
    {
        _params.insert({"edge", edge}); 
        this->param_rangecheck("edge"); 
    }
    else if (name == "cuboid")
    {
        vector<double> my_edges (_dim, 1.0);
        if (edges.size() != 0 && edges.size() != _dim)
            {throw std::invalid_argument("Cuboid's |edges| neither 0 nor _dim");}
        else if (edges.size() == 0)
            {edges = my_edges;}
        for (int i = 0; i<_dim; i++)
        {
            double edge_i = edges[i];
            std::string key_s = "edge_"+std::to_string(i);
            const char* key_i = key_s.c_str();
            //char key_i[key_s.length() + 1]; 
	        //std::strcpy(key_i, key_s.c_str());
            _params.insert({key_i, edge_i}); 
            this->param_rangecheck(key_i);
        } 
    }
}
       
        
void CoordinateShape::param_rangecheck(const char* name, 
                                        double maxValue,
                                        bool canBeZero)
{
    double p = _params.find(name)->second;
    bool isTooLow = ((canBeZero && p < 0.0) or (!canBeZero && p<=0));
    const char* errorchar = (canBeZero) ? "greater than or equal to 0" 
                                        : "greater than 0";
    if (std::isnan(maxValue) and isTooLow)
    {
        std::string error = "Parameter ";
        error += name;
        error += " is out of range. It must be a double ";
        error += errorchar;
        throw std::invalid_argument(error);
    }
    else if (isTooLow or p>=maxValue)
    {
        std::string error = "Parameter ";
        error += name;
        error += " is out of range. It must be a double ";
        error += errorchar;
        error += " and smaller than ";
        error += std::to_string(maxValue);
        throw std::invalid_argument(error);
    }

}


double CoordinateShape::Parameter (const char* key)
/**
 * @brief Return value of parameter corresponding to "key". 
 *        Note: for cuboid, need to call "edge_i" where i in [0, _dim-1].
 */
    {return _params.find(key)->second;}


vector<double> CoordinateShape::Edges ()
/**
 * @brief Return value of [0, _dim-1] edges.
 * @exception Raise error if undefined.
 */
{
    vector<double> edges (_dim);
    for (int i = 0; i<_dim; i++)
    {
        std::string key_s = "edge_"+std::to_string(i);
        const char* key_i = key_s.c_str();
        edges[i] = _params.find(key_i)->second;
    }
    return edges;
}


double CoordinateShape::Volume()
/**
 * @brief Return Volume. On first call volume is computed. 
 */
{
    if (_volume != 0)
        {return _volume;}
    else if (_name == "ball")
    {
        double r = this->Parameter("radius");
        _volume  = std::pow(r, _dim); //for now
        double h_fr = this->Parameter("hollow");
        _volume -= (h_fr>0)? std::pow(h_fr*r, _dim) : 0.0;
        _volume *= std::pow(M_PI, _dim/2) / std::tgamma(_dim/2 +1);
    }
    else if (_name == "cylinder")
    {
        double r = this->Parameter("radius");
        _volume  = std::pow(r, _dim-1); //for now
        double h_fr = this->Parameter("hollow");
        _volume -= (h_fr>0)? std::pow(h_fr*r, _dim-1) : 0.0;
        _volume *= std::pow(M_PI, _dim/2 -0.5) / std::tgamma(_dim/2 +0.5);
        _volume *= this->Parameter("duration");
    }
    else if (_name == "bicone")
    {
        double r = this->Parameter("radius");
        _volume  = std::pow(r, _dim); //for now
        double h_fr = this->Parameter("hollow");
        _volume -= (h_fr>0)? std::pow(h_fr*r, _dim) : 0.0;
        _volume *= std::pow(M_PI, _dim/2 -0.5) / std::tgamma(_dim/2 +0.5);
        _volume *= 2/_dim;
    }
    else if (_name == "cube")
    {
        _volume = std::pow(this->Parameter("edge"), _dim);
    }
    else if (_name == "cuboid")
    {
        _volume = 1;
        for (int i = 0; i<_dim; i++)
        {
            std::string key_s = "edge_"+std::to_string(i);
            const char* key_i = key_s.c_str();
            _volume *= this->Parameter(key_i);
        }
    }
    return _volume;
}


vector<double> CoordinateShape::Limits(int dim)
/**
 * @brief Returns a vector<double> for the minimum (left) and maximum (right) 
 *        of a shape along coordinate dimension "dim" (0-indexed).
 */
{
    if ((dim < 0) || (dim >= _dim))
    {
        throw std::invalid_argument("Dimension " +std::to_string(dim) +
                                    "out of range " + "[0, " 
                                    + std::to_string(_dim)+" ]!");
    }
    double l = 0;
    if ((dim == 0) && (_name == "cylinder"))
        {l = this->Parameter("duration") /2;}
    else if (_name == "cube")
        {l = this->Parameter("edge") /2;}
    else if (_name == "cuboid")
        {l = this->Parameter(("edge_"+std::to_string(dim)).c_str()) /2;}
    else
        {l = this->Parameter("radius");}
    double shift = _center[dim];
    vector<double> limits = {-l + shift, l + shift};
    return limits; 
}

//CoordinateShape::~CoordinateShape() {}

// int main(){

//     std::cout << "shapes.cpp WORKS!" << std::endl;
// }