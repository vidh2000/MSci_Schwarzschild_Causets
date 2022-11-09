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
#include <string.h>
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
       {std:: cout << "Dim smaller than 1!" << std::endl;
        throw std::invalid_argument("Dim smaller than 1!");}
    _dim = dim;

    // Set Name
    isBicone = strcmp(name,"diamond")==0 || strcmp(name, "bicone")==0;
    if (strcmp(name,"ball")!=0 && !isBicone &&
        strcmp(name,"cylinder")!=0 && strcmp(name,"cube")!=0 &&
        strcmp(name,"cuboid")!=0)
     { 
        std::cout << "The given shape is "<< name << ". Choose: ";
        std::cout << "'ball', 'bicone', 'cylinder', 'cube' or 'cuboid'.";
        std::cout << std::endl; 
        throw std::invalid_argument("Wrong name chosen");
    }
     if (isBicone){
        _name = "bicone";}
     else{
        _name = name;}

    // Set Center
    vector<double> zero_center(_dim, 0.0);
    if (center.size() == dim)
        {_center = center;}
    else if (center.size() == 0){
        _center = zero_center;}
    else {
        std::cout << "Center's size neither 0 nor _dim" << std::endl;
        throw std::invalid_argument("Center's size neither 0 nor _dim");}
    
    std::cout << "isBicone before radius setting = " << isBicone << std::endl;

    // Set Shape Parameters
    if (strcmp(name,"ball")==0 || isBicone)
    {
        std::cout << "radius = " << radius << std::endl;
        _params["radius"] = radius;//_params.insert({"radius", radius}); 
        this->param_rangecheck("radius");
        std::cout << "radius param = " << _params["radius"] << std::endl;

        _params["hollow"] = hollow;//.insert({"hollow", hollow}); 
        this->param_rangecheck("hollow", 1.0, true);
        std::cout << "radius param after hollow = "<< _params["radius"] << std::endl;
    }
    else if (strcmp(name,"cylinder")==0)
    {
        _params.insert({"radius", radius}); 
        this->param_rangecheck("radius");

        _params.insert({"hollow", hollow}); 
        this->param_rangecheck("hollow", 1.0, true); 

        _params.insert({"duration", duration}); 
        this->param_rangecheck("duration");   
    }
    else if (strcmp(name,"cube")==0)
    {
        _params.insert({"edge", edge}); 
        this->param_rangecheck("edge"); 
    }
    else if (strcmp(name,"cuboid")==0)
    {
        vector<double> my_edges (_dim, 1.0);
        if (edges.size() != 0 && edges.size() != _dim){
            std::cout << "Cuboid's |edges| neither 0 nor _dim" << std::endl;
            throw std::invalid_argument("Cuboid's |edges| neither 0 nor _dim");}
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
    std::cout << "Inside param_rangecheck: name = " << name << std::endl;
    bool isTooLow = ((canBeZero && p < 0.0) or (!canBeZero && p<=0));
    const char* errorchar = (canBeZero) ? "greater than or equal to 0" 
                                        : "greater than 0";
    if (std::isnan(maxValue) and isTooLow)
    {
        std::string error = "Parameter ";
        error += name;
        error += " is out of range. It must be a double ";
        error += errorchar;
        std::cout << error << std::endl;
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
        std::cout << error << std::endl;
        throw std::invalid_argument(error);
    }

}


double CoordinateShape::Parameter (const char* key)
/**
 * @brief Return value of parameter corresponding to "key". 
 *        Note: for cuboid, need to call "edge_i" where i in [0, _dim-1].
 */
    {return _params[key];} 


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
    std::cout << "in CoordinateShape::Volume(), isBicone =" << isBicone << std::endl;
    if (_volume != 0)
        {return _volume;}
    else if (strcmp(_name, "ball")==0)
    {
        double r = this->Parameter("radius");
        _volume  = std::pow(r, _dim); //for now
        double h_fr = this->Parameter("hollow");
        _volume -= (h_fr>0)? std::pow(h_fr*r, _dim) : 0.0;
        _volume *= std::pow(M_PI, _dim/2) / std::tgamma(_dim/2 +1);
    }
    else if (strcmp(_name, "cylinder")==0)
    {
        double r = this->Parameter("radius");
        _volume  = std::pow(r, _dim-1); //for now
        double h_fr = this->Parameter("hollow");
        _volume -= (h_fr>0)? std::pow(h_fr*r, _dim-1) : 0.0;
        _volume *= std::pow(M_PI, _dim/2 -0.5) / std::tgamma(_dim/2 +0.5);
        _volume *= this->Parameter("duration");
    }
    else if (isBicone)
    {
        double r = this->Parameter("radius");
        std::cout << "r = " << r << std::endl;
        _volume  = std::pow(r, _dim); //for now
        double h_fr = this->Parameter("hollow");
        _volume -= (h_fr>0)? std::pow(h_fr*r, _dim) : 0.0;
        _volume *= std::pow(M_PI, _dim/2 -0.5) / std::tgamma(_dim/2 +0.5);
        _volume *= 2/_dim;
    }
    else if (strcmp(_name, "cube")==0)
    {
        _volume = std::pow(this->Parameter("edge"), _dim);
    }
    else if (strcmp(_name,"cuboid")==0)
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
        std::string msg = "Dimension " +std::to_string(dim) +
                                    "out of range " + "[0, " 
                                    + std::to_string(_dim)+" ]!";
        std::cout << msg << std::endl;
        throw std::invalid_argument("Dimension " +std::to_string(dim) +
                                    "out of range " + "[0, " 
                                    + std::to_string(_dim)+" ]!");
    }
    double l = 0;
    if ((dim == 0) && (strcmp(_name, "cylinder")==0))
        {l = this->Parameter("duration") /2;}
    else if (strcmp(_name,"cube")==0)
        {l = this->Parameter("edge") /2;}
    else if (strcmp(_name, "cuboid")==0)
        {l = this->Parameter(("edge_"+std::to_string(dim)).c_str()) /2;}
    else
        {l = this->Parameter("radius");}
    double shift = _center[dim];
    vector<double> limits = {-l + shift, l + shift};
    return limits; 
}
//CoordinateShape::~CoordinateShape() {}


// int main(){
//     CoordinateShape S;
//     std::cout << "S._name = " << S._name << std::endl;
//     std::cout << "radius = " << S.Parameter("radius") << std::endl;
//     std::cout << "shapes.cpp WORKS!" << std::endl;
// }