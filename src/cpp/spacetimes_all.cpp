#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <stdio.h>
#include <string.h>
#include <stack>
#include <string>
#include <vector>

#include "boost/numeric/odeint.hpp"

#include "shapes.h"
#include "spacetimes.h"
#include "vecfunctions.h"
#include "functions.h"

using std::vector;

#define M_PI  3.14159265358979323846  /* pi */


/**
 * @brief Construct a new Spacetime:: Spacetime object. It does nothing 
 * but initialise a Spacetime onject. Then call:
 * - FlatSpacetime(dim, period)
 * - BlackHoleSpacetime(dim)
 */
Spacetime::Spacetime(){}


/**
 * @brief Internal function for the time sampling array for a cone from 
 * "origin point"(hence take t0=origin[0]) to time "t".
 */
vector<double> Spacetime::T_slice_sampling(double t, 
                                            vector<double>origin,
                                            int samplingsize)
{
    vector<double> ts (samplingsize);
    double t_i = (t-origin[0])/samplingsize;
    for (int i = 0; i<samplingsize; i++)
        {ts[i] += origin[0] + t_i*i;}
    return ts;
}     



typedef bool (*func)
(const vector<double>& xvec, const vector<double>& yvec, 
 vector<double> period, double mass);
/**
 * @brief Return callable "bool callable (const vector<double>& xvec, 
 * const vector<double>& yvec, vector<double> period, double mass", i.e.
 * returns a function that returns the bool saying if x-y is timelike
 * based on spacetime features.
 */
func Spacetime::Causality()
{
    if (_dim == 1)
        {return &Spacetime::causal1d;}
    else if (std::strcmp(_name, "Flat")==0)
    {
        if (_isPeriodic) return &Spacetime::Flat_causal_periodic;
        else return &Spacetime::Flat_causal;
    }
    else /*(std::strcmp(_name, "BlackHole")==0)*/
    {
        // Changing coordinates into EF (original)
        // place holder //
        if      (_dim == 4) return &Spacetime::BH_causal4D;
        else if (_dim == 3) return &Spacetime::BH_causal3D;
        else if (_dim == 2) return &Spacetime::BH_causal2D;
        else
        {
            std::cout << "Please use D=2, D=3 or D=4 for Schwarzschild\n";
            throw std::invalid_argument("problem");
        }
    }
}


typedef vector<bool> (*general_func)
(const vector<double>& xvec, const vector<double>& yvec, vector<double> period, double mass);
/**
 * @brief Return callable "vector<bool> callable (const vector<double>& xvec, 
 * const vector<double>& yvec, vector<double> period, double mass", i.e.
 * returns a function that returns a vector of bools {x-y timelike, x<=y, x>y},
 * based on spacetime features.
 */
general_func Spacetime::General_Causality()
{
    if (_dim == 1)
        {return &Spacetime::general_causal1d;}
    else if (std::strcmp(_name, "Flat")==0)
    {
        if (_isPeriodic) return &Spacetime::Flat_general_causal_periodic;
        else return &Spacetime::Flat_general_causal;
    }
    else /*(std::strcmp(_name, "BlackHole")==0)*/
    {
        // Changing coordinates into EF (original)
        // place holder //
        if      (_dim == 4) return &Spacetime::BH_general_causal4D;
        else if (_dim == 3) return &Spacetime::BH_general_causal3D;
        else if (_dim == 2) return &Spacetime::BH_general_causal2D;
        else
        {
            std::cout << "Please use D=2, D=3 or D=4 for Schwarzschild\n";
            throw std::invalid_argument("problem");
        }
    }
}


/**
 * @brief Simplest function of two events in 1D
 * 
 * @param xvec : vector<double>(1), time-coordinate of x
 * @param yvec : vector<double>(1), time-coordinate of y
 * @param period : needed for consistency with Causality, not used
 * @return bool : 1
 */
bool Spacetime::causal1d(const vector<double>& xvec, const vector<double>& yvec,
                         vector<double> period, double mass)
{
    return true;
}


/**
 * @brief Simplest function of two events in 1D
 * 
 * @param xvec : vector<double>(1), time-coordinate of x
 * @param yvec : vector<double>(1), time-coordinate of y
 * @param period : needed for consistency with Causality, not used
 * @return vector<bool> : {1, x<=y, x>y}
 */
vector<bool> Spacetime::general_causal1d(const vector<double>& xvec, 
                                         const vector<double>& yvec,
                                         vector<double> period, double mass)
{
    //std::cout<<"For debuggin: calling Spacetime::causal1d\n";
    double t_delta = yvec[0] - xvec[0];
    return {1, t_delta >= 0.0, t_delta < 0.0};
}


/**
 * @brief Convert from cartesian to spherical coordinates.
 */
void Spacetime::CarttoSpherical (std::vector<double>& xvec)
{
    if (xvec.size()==3)
    {
        double r = std::sqrt(xvec[1]*xvec[1] + xvec[2]*xvec[2]);
        double phi = std::atan2(xvec[2],xvec[1]);
        xvec[1] = r;
        xvec[2] = phi;
    }
    else if (xvec.size()==4)
    {
        double rho = std::sqrt(xvec[1]*xvec[1] + xvec[2]*xvec[2]);
        double r = std::sqrt(rho*rho + xvec[3]*xvec[3]);
        double theta = std::atan2(rho, xvec[3]);
        double phi = std::atan2(xvec[2],xvec[1]);
        xvec[1] = r;
        xvec[2] = theta;
        xvec[3] = phi;
    }
}


/**
 * @brief Convert from cartesian to spherical coordinates.
 */
void Spacetime::CarttoSpherical (std::vector<std::vector<double>>& coords)
{
    if (coords[0].size()==3)
    {
        for (auto & xvec : coords)
        {
            double r = std::sqrt(xvec[1]*xvec[1] + xvec[2]*xvec[2]);
            double phi = std::atan2(xvec[2],xvec[1]);
            xvec[1] = r;
            xvec[2] = (phi>0)? phi : phi + 2*M_PI;
        }
    }
    else if (coords[0].size()==4)
    {
        for (auto & xvec : coords)
        {
            double rho = std::sqrt(xvec[1]*xvec[1] + xvec[2]*xvec[2]);
            double r = std::sqrt(rho*rho + xvec[3]*xvec[3]);
            double theta = std::atan2(rho, xvec[3]);
            double phi = std::atan2(xvec[2],xvec[1]);
            xvec[1] = r;
            xvec[2] = theta;
            xvec[3] = (phi>0)? phi : phi + 2*M_PI;
        }
    }
}

/**
 * @brief Convert spherical to cartesian coordinates
 * 
 * @param coords coordinates in spherical coord.syst (t,r,theta,phi) 
 */
void Spacetime::SphericaltoCart(std::vector<std::vector<double>>& coords)
{
    if (coords[0].size()==3)
    {
        for (auto & xvec : coords)
        {
            double x = xvec[1]*std::cos(xvec[2]);
            double y = xvec[1]*std::sin(xvec[2]);
            xvec[1] = x;
            xvec[2] = y;
        }
    }
    else if (coords[0].size()==4)
    {
        for (auto & xvec : coords)
        {
            double x = xvec[1]*std::sin(xvec[2])*std::cos(xvec[3]);
            double y = xvec[1]*std::sin(xvec[2])*std::sin(xvec[3]);
            double z = xvec[1]*std::cos(xvec[2]);
            xvec[1] = x;
            xvec[2] = y;
            xvec[3] = z;
        }
    }
}


/**
 * @brief Turn ingoing Eddington Finkelstein coordinates into Schwarzschild.
 * 
 * @param xvec std::vector<double>& : vector which gets changed
 * @param mass double : mass of BH
 * @param EFtype const char* :
 * - "original" : t_EF = ts + 2M ln|r/2M - 1| (as in He, Ridoeut)
 * - "uv" : t_EF = ts + r + 2M ln|r/2M - 1| (as uv in Wikipedia's EF)
 */
void Spacetime::InEFtoS (std::vector<double> &xvec, double mass,
                        const char* EFtype)
{
    if (strcmp(EFtype, "original")==0)
    {
        double arg = std::abs(xvec[1])/(2*mass) - 1;
        if (arg>=0) xvec[0] -= 2*mass*std::log(arg);
        else        xvec[0] -= 2*mass*std::log(-arg);
    }
    else if (strcmp(EFtype, "uv")==0)
    {
        double arg = std::abs(xvec[1])/(2*mass) - 1;
        if (arg>=0) xvec[0] -= std::abs(xvec[1]) + 2*mass*std::log(arg);
        else        xvec[0] -= std::abs(xvec[1]) + 2*mass*std::log(-arg);
    }
    else
    {
        std::cout<<"method in StoinEF must be 'original' or 'uv'\n";
        throw std::invalid_argument("Wrong method");
    }
}

/**
 * @brief Turn ingoing Eddington Finkelstein coordinates into Schwarzschild.
 * 
 * @param coords std::vector<vector<double>>& : list of coordinates getting
 *  changed
 * @param mass double : mass of BH
 * @param EFtype const char* :
 * - "original" : t_EF = ts + 2M ln|r/2M - 1| (as in He, Ridoeut)
 * - "uv" : t_EF = ts + r + 2M ln|r/2M - 1| (as uv in Wikipedia's EF)
 */
void Spacetime::InEFtoS (std::vector<std::vector<double>> &coords, double mass,
                        const char* EFtype)
{
    for (std::vector<double> & xvec : coords)
        {Spacetime::InEFtoS(xvec, mass, EFtype);}
}


/**
 * @brief Turn Schwarzschild coordinates into ingoing Eddington Finkelstein.
 * 
 * @param xvec std::vector<double>& : vector getting changed
 * @param mass double : mass of BH
 * @param EFtype const char* :
 * - "original" : t_EF = ts + 2M*ln|r/2M - 1| (as in He, Ridoeut)
 * - "uv" : t_EF = ts + r + 2M*ln|r/2M - 1| (as uv in Wikipedia's EF)
 */
void Spacetime::StoInEF (std::vector<double> &xvec, double mass,
                        const char* EFtype)
{
    if (strcmp(EFtype, "original")==0)
    {
        double arg = std::abs(xvec[1])/(2*mass) - 1;
        if (arg>=0) xvec[0] += 2*mass*std::log(arg);
        else        xvec[0] += 2*mass*std::log(-arg);
    }
    else if (strcmp(EFtype, "uv")==0)
    {
        double arg = std::abs(xvec[1])/(2*mass) - 1;
        if (arg>=0) xvec[0] += std::abs(xvec[1]) + 2*mass*std::log(arg);
        else        xvec[0] += std::abs(xvec[1]) + 2*mass*std::log(-arg);
    }
    else
    {
        //std::cout<<"method in StoinEF must be 'original' or 'uv'\n";
        throw std::invalid_argument("Wrong method");
    }
}


/**
 * @brief Turn Schwarzschild coordinates into ingoing Eddington Finkelstein.
 * 
 * @param coords vector<vector<double>>& : list of coordinates getting changed
 * @param mass double : mass of BH
 * @param EFtype const char* :
 * - "original" : t_EF = ts + 2M*ln|r/2M - 1| (as in He, Ridoeut)
 * - "uv" : t_EF = ts + r + 2M*ln|r/2M - 1| (as uv in Wikipedia's EF)
 */
void Spacetime::StoInEF (std::vector<std::vector<double>> &coords, double mass,
                        const char* EFtype)
{
    for (std::vector<double> & xvec : coords)
        {Spacetime::StoInEF(xvec, mass, EFtype);}
}


/**
 * @brief Turn Gullstrand–Painlevé into Schwarzschild coordinates.
 * 
 * @param coords vector<double>& : vector getting changed
 * @param mass double : mass of BH
 */
void Spacetime::GPtoS (std::vector<double>& xvec,  double mass)
{
    double y = std::sqrt(xvec[1]/(2*mass));
    xvec[0] -= 2*mass*( -2*y + std::log( (y+1)/(y-1)) );
}


/**
 * @brief Turn Gullstrand–Painlevé into Schwarzschild coordinates.
 * 
 * @param coords vector<vector<double>>& : list of coordinates getting changed
 * @param mass double : mass of BH
 */
void Spacetime::GPtoS (std::vector<std::vector<double>>& coords, 
                       double mass)
{
    for (std::vector<double> & xvec : coords)
        {Spacetime::GPtoS(xvec, mass);}
}


/**
 * @brief Turn Schwarzschild coordinates into Gullstrand–Painlevé.
 * 
 * @param coords vector<double>& : vector getting changed
 * @param mass double : mass of BH
 */
void Spacetime::StoGP (std::vector<double>& xvec, double mass)
{
    double y = std::sqrt(xvec[1]/(2*mass));
    xvec[0] += 2*mass*( -2*y + std::log( (y+1)/(y-1)) );
}


/**
 * @brief Turn Schwarzschild coordinates into Gullstrand–Painlevé.
 * 
 * @param coords vector<vector<double>>& : list of coordinates getting changed
 * @param mass double : mass of BH
 */
void Spacetime::StoGP (std::vector<std::vector<double>>& coords, 
                       double mass)
{
    for (std::vector<double> & xvec : coords)
        {Spacetime::StoGP(xvec, mass);}
}


/**
 * @brief Turn ingoing Eddington Finkelstein coordinates into 
 * Gullstrand–Painlevé.
 * 
 * @param xvec std::vector<double>& : vector which gets changed
 * @param mass double : mass of BH
 * @param EFtype const char* :
 * - "original" : t_EF = ts + 2M ln|r/2M - 1| (as in He, Ridoeut)
 * - "uv" : t_EF = ts + r + 2M ln|r/2M - 1| (as uv in Wikipedia's EF)
 */
void Spacetime::InEFtoGP (std::vector<double>& xvec, double mass,
                          const char* EFtype)
{
    if (strcmp(EFtype, "original")==0)
    {
        double y = std::sqrt(xvec[1]/(2*mass));
        xvec[0] -= 2*mass*( -2*y + std::log( (y+1)/(y-1)) );
        double arg = std::abs(xvec[1])/(2*mass) - 1;
        if (arg>=0) xvec[0] -= 2*mass*std::log(arg);
        else        xvec[0] -= 2*mass*std::log(-arg);
    }
    else if (strcmp(EFtype, "uv")==0)
    {
        double y = std::sqrt(xvec[1]/(2*mass));
        xvec[0] -= 2*mass*( -2*y + std::log( (y+1)/(y-1)) );
        double arg = std::abs(xvec[1])/(2*mass) - 1;
        if (arg>=0) xvec[0] -= std::abs(xvec[1]) + 2*mass*std::log(arg);
        else        xvec[0] -= std::abs(xvec[1]) + 2*mass*std::log(-arg);
    }
    else
    {
        //std::cout<<"method in StoinEF must be 'original' or 'uv'\n";
        throw std::invalid_argument("Wrong method");
    }
}


/**
 * @brief Turn ingoing Eddington Finkelstein coordinates into 
 * Gullstrand–Painlevé.
 * 
 * @param xvec std::vector<vector<double>>& : list of vectors geting changed
 * @param mass double : mass of BH
 * @param EFtype const char* :
 * - "original" : t_EF = ts + 2M ln|r/2M - 1| (as in He, Ridoeut)
 * - "uv" : t_EF = ts + r + 2M ln|r/2M - 1| (as uv in Wikipedia's EF)
 */
void Spacetime::InEFtoGP (std::vector<std::vector<double>>& coords, 
                          double mass, const char* EFtype)
{
    for (std::vector<double> & xvec : coords)
        {InEFtoGP(xvec, mass, EFtype);}
}


/**
 * @brief Turn Gullstrand–Painlevé into ingoing Eddington Finkelstein 
 * coordinates.
 * 
 * @param xvec std::vector<double>& : vector which gets changed
 * @param mass double : mass of BH
 * @param EFtype const char* :
 * - "original" : t_EF = ts + 2M ln|r/2M - 1| (as in He, Ridoeut)
 *-----------------------------------------------------------------------------
 * - "uv" : t_EF = ts + r + 2M ln|r/2M - 1| (as uv in Wikipedia's EF)
 */
void Spacetime::GPtoInEF (std::vector<double>& xvec, double mass,
                          const char* EFtype)
{
    if (strcmp(EFtype, "original")==0)
    {
        double y = std::sqrt(xvec[1]/(2*mass));
        xvec[0] += 2*mass*( -2*y + std::log( (y+1)/(y-1)) );
        double arg = std::abs(xvec[1])/(2*mass) - 1;
        if (arg>=0) xvec[0] += 2*mass*std::log(arg);
        else        xvec[0] += 2*mass*std::log(-arg);
    }
    else if (strcmp(EFtype, "uv")==0)
    {
        double y = std::sqrt(xvec[1]/(2*mass));
        xvec[0] += 2*mass*( -2*y + std::log( (y+1)/(y-1)) );
        double arg = std::abs(xvec[1])/(2*mass) - 1;
        if (arg>=0) xvec[0] += std::abs(xvec[1]) + 2*mass*std::log(arg);
        else        xvec[0] += std::abs(xvec[1]) + 2*mass*std::log(-arg);
    }
    else
    {
        //std::cout<<"method in StoinEF must be 'original' or 'uv'\n";
        throw std::invalid_argument("Wrong method");
    }
}


/**
 * @brief Turn Gullstrand–Painlevé into ingoing Eddington Finkelstein 
 * coordinates.
 * 
 * @param xvec std::vector<vector<double>>& : list of vectors geting changed
 * @param mass double : mass of BH
 * @param EFtype const char* :
 * - "original" : t_EF = ts + 2M ln|r/2M - 1| (as in He, Ridoeut)
 *-----------------------------------------------------------------------------
 * - "uv" : t_EF = ts + r + 2M ln|r/2M - 1| (as uv in Wikipedia's EF)
 */
void Spacetime::GPtoInEF (std::vector<std::vector<double>>& coords,
                          double mass, const char* EFtype)
{
    for (std::vector<double> & xvec : coords)
        {GPtoInEF(xvec, mass, EFtype);}
}


/**
 * @brief Turn ingoing Eddington Finkelstein coordinates into KS.
 * 
 * @param xvec std::vector<double>& : vector which gets changed
 * @param mass double : mass of BH
 * @param EFtype const char* :
 * - "original" : t_EF = ts + 2M ln|r/2M - 1| (as in He, Ridoeut)
 *-----------------------------------------------------------------------------
 * - "uv" : t_EF = ts + r + 2M ln|r/2M - 1| (as uv in Wikipedia's EF)
 */
void Spacetime::InEFtoKS (std::vector<double>& xvec, double mass,
                          const char* EFtype)
{
    // Define this for later
    double arg = std::abs(xvec[1])/(2*mass) - 1;

    // First turn into S
    if (strcmp(EFtype, "original")==0)
    {
        if (arg>=0) xvec[0] -= 2*mass*std::log(arg);
        else        xvec[0] -= 2*mass*std::log(-arg);
    }
    else if (strcmp(EFtype, "uv")==0)
    {
        if (arg>=0) xvec[0] -= std::abs(xvec[1]) + 2*mass*std::log(arg);
        else        xvec[0] -= std::abs(xvec[1]) + 2*mass*std::log(-arg);
    }
    else
    {
        std::cout<<"method in StoinEF must be 'original' or 'uv'\n";
        throw std::invalid_argument("Wrong method");
    }

    // Then turn into X and T
    double T = std::sqrt(arg) * std::exp(xvec[1]/(4*mass))
            *std::sinh(xvec[0]/(4*mass));
    double X = std::sqrt(arg) * std::exp(xvec[1]/(4*mass))
            *std::cosh(xvec[0]/(4*mass));
    xvec[0] = T;
    xvec[1] = X;
}


/**
 * @brief Turn ingoing Eddington Finkelstein coordinates into KS.
 * 
 * @param xvec std::vector<vector<double>>& : list of vectors geting changed
 * @param mass double : mass of BH
 * @param EFtype const char* :
 * --- "original" : t_EF = ts + 2M ln|r/2M - 1| (as in He, Ridoeut
 *-----------------------------------------------------------------------------)
 * --- "uv" : t_EF = ts + r + 2M ln|r/2M - 1| (as uv in Wikipedia's EF)
 */
void Spacetime::InEFtoKS (std::vector<std::vector<double>>& coords, 
                          double mass, const char* EFtype)
{
    for (std::vector<double> & xvec : coords)
        {InEFtoKS(xvec, mass, EFtype);}
}


/**
 * @brief Turn KS into ingoing Eddington Finkelstein coordinates.
 * 
 * @param xvec std::vector<double>& : vector which gets changed
 * @param mass double : mass of BH
 * @param EFtype const char* :
 * - "original" : t_EF = ts + 2M ln|r/2M - 1| (as in He, Ridoeut)
 *-----------------------------------------------------------------------------
 * - "uv" : t_EF = ts + r + 2M ln|r/2M - 1| (as uv in Wikipedia's EF)
 */
void Spacetime::KStoInEF (std::vector<double>& xvec, double mass,
                          const char* EFtype)
{
    //Turn from KS to S
    double t = 0;
    double r = 0;
    double delta2 = std::pow(xvec[0],2) - std::pow(xvec[1],2);
    if (delta2 < 0 && xvec[1] > 0)
    {
        std::cout<<"I am sorry, currently not implemented";
        throw std::invalid_argument("Not impplemented");
    }
    else if (1 > delta2 && delta2 > 0 && xvec[0] > 0)
    {
        std::cout<<"I am sorry, currently not implemented";
        throw std::invalid_argument("Not impplemented");
    }

    // Turn from S to inEF
    if (strcmp(EFtype, "original")==0)
    {
        double arg = std::abs(xvec[1])/(2*mass) - 1;
        if (arg>=0) xvec[0] += 2*mass*std::log(arg);
        else        xvec[0] += 2*mass*std::log(-arg);
    }
    else if (strcmp(EFtype, "uv")==0)
    {
        double arg = std::abs(xvec[1])/(2*mass) - 1;
        if (arg>=0) xvec[0] += std::abs(xvec[1]) + 2*mass*std::log(arg);
        else        xvec[0] += std::abs(xvec[1]) + 2*mass*std::log(-arg);
    }
    else
    {
        std::cout<<"method in StoinEF must be 'original' or 'uv'\n";
        throw std::invalid_argument("Wrong method");
    }
}


/**
 * @brief Turn KS into ingoing Eddington Finkelstein coordinates.
 * 
 * @param xvec std::vector<vector<double>>& : list of vectors geting changed
 * @param mass double : mass of BH
 * @param EFtype const char* :
 * --- "original" : t_EF = ts + 2M ln|r/2M - 1| (as in He, Ridoeut)
 * ---------------------------------------------------------------------------
 * --- "uv" : t_EF = ts + r + 2M ln|r/2M - 1| (as uv in Wikipedia's EF)
 */
void Spacetime::KStoInEF (std::vector<std::vector<double>>& coords,
                          double mass, const char* EFtype)
{
    for (std::vector<double> & xvec : coords)
        {KStoInEF(xvec, mass, EFtype);}
}



/**
 * @brief Switch type of EF coordinate from/to original to/from uv
 * 
 * @param coords vector<vector<double>> : one vector of coordinates to change
 * @param from const char* : from which type of coordinate to change
 * - "original" : then goes to "uv";
 *-----------------------------------------------------------------------------
 * - "uv" : then goes to "original";
 */
void switchInEF (std::vector<double>& xvec, const char* from)
{
    if (strcmp(from, "original")==0) xvec[0] += std::abs(xvec[1]);
    else xvec[0] -= std::abs(xvec[1]);
}


/**
 * @brief Switch type of EF coordinate
 * 
 * @param coords vector<vector<double>> : list of coordinates to change
 * @param from const char* : from which type of coordinate to change
 * - "original" : then goes to "uv";
 *-----------------------------------------------------------------------------
 * - "uv" : then goes to "original";
 */
void Spacetime::switchInEF (std::vector<std::vector<double>>& coords, 
                            const char* from)
{
    if (strcmp(from, "original")==0)
    {
        for (std::vector<double> & xvec : coords)
            {xvec[0] += std::abs(xvec[1]);}
    }
    else
    {
        for (std::vector<double> & xvec : coords)
            {xvec[0] -= std::abs(xvec[1]);}
    }
}




///////////////////////////////////////////////////////////////////////////////
/// GENERAL CAUSALITIES ///////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/**
 * @brief Causality algorithm for two events in 2D EForig coordinates, from
 * Song He and David Rideout 2009 Class. Quantum Grav. 26 125015. 
 * 
 * @param xvec vector<double> : EF coordinates of x.
 * @param yvec vector<double> : EF coordinates of y.
 * @param period vector<double> : period along SPATIAL coordinates. Currently 
 * not implemented in BH, hence deafult is {}.
 * @param mass : mass of Black Hole
 * @return vector<bool> : {x-y timelike, x<=y, x>y}
 */
vector<bool> Spacetime::BH_general_causal2D (const vector<double>& xvec, 
                                            const vector<double>& yvec,
                                            std::vector<double> period,
                                            double mass)
{
    //WORKING IN EF ORIGINAL COORDINATES
    if (yvec[0]<xvec[0])
    {
        std:vector<bool> result = Spacetime::BH_general_causal2D
                                             (yvec, xvec, period);
        if (result[0])
        {
            bool a = result[1]*1;
            result[1] = result[2]*1;
            result[2] = result[1]*1;
        }
        return result;
    }

    double t1     = xvec[0]; double t2     = yvec[0];
    double r1     = xvec[1]; double r2     = yvec[1];

    if (r1*r2 < 0) //on opposite sides of the singularity
        {return {false, false, false};}
    else
    {
        r1 = std::abs(r1);
        r2 = std::abs(r2);
    }
    
    // Section 2.2: Radially separated pairs and radial null geodesics
    // as in 2D necessarily varphi2 = 0
    if (r1>=r2)
    {
        // all 3 cases of the paper return same
        if (t2 >= t1 + r1 - r2 - 1e-6)
            {return {true, true, false};}
        else
            {return {false, false, false};}
    }
    else /*r2>r1*/
    {
        if (r1<=2*mass)
            {return {false, false, false};}
        else // if (r1>2*mass)
        {
            if (t2 >= t1 + r2 - r1 + 4*mass*std::log((r2-2*mass)/(r1-2*mass)))
                {return {true, true, false};}
            else
                {return {false, false, false};}
        }
    }
}


/**
 * @brief Causality algorithm for two events in 4D EForig coordinates, from
 * Song He and David Rideout 2009 Class. Quantum Grav. 26 125015. 
 * 
 * @param xvec vector<double> : EF coordinates of x.
 * @param yvec vector<double> : EF coordinates of y.
 * @param period vector<double> : period along SPATIAL coordinates. Currently 
 * not implemented in BH, hence deafult is {}.
 * @param mass : mass of Black Hole
 * @return vector<bool> : {x-y timelike, x<=y, x>y}
 */
vector<bool> Spacetime::BH_general_causal3D (const vector<double>& xvec, 
                                    const vector<double>& yvec,
                                    std::vector<double> period,
                                    double mass)
{
    //IF WORKING IN EF COORDINATES
    if (yvec[0]<xvec[0])
    {
        std:vector<bool> result = Spacetime::BH_general_causal3D
                                            (yvec, xvec, period);
        if (result[0])
        {
            bool a = result[1]*1;
            result[1] = result[2]*1;
            result[2] = result[1]*1;
        }
        return result;
    }

    double t1     = xvec[0]; double t2     = yvec[0];
    double r1     = xvec[1]; double r2     = yvec[1];
    double phi1   = xvec[2]; double phi2   = yvec[2];

    double vartheta1 = M_PI / 2;
    double vartheta2 = M_PI / 2; 
    double varphi1 = 0;
    double varphi2 = phi1-phi2;
    if (varphi2 < 0)
        {varphi2 += 2*M_PI;}
    
    vector<double> transf_xvec = {t1, r1, vartheta1, varphi1};
    vector<double> transf_yvec = {t2, r2, vartheta2, varphi2};
    
    // Section 2.2: Radially separated pairs and radial null geodesics
    if (varphi2<1e-6) //should be ==zero, but leave room for some error
    {
        if (r1>=r2)
        {
            // all 3 cases of the paper return same
            if (t2 >= t1 + r1 - r2 - 1e-6)
                    {return {true, true, false};}
            else
                {return {false, false, false};}
        }
        else
        {
            if (r1<=2*mass)
                {return {false, false, false};}
            else // if (r1>2*mass)
            {
                if (t2 >= t1 + r2 - r1 + 
                            4*mass*std::log((r2-2*mass)/(r1-2*mass)))
                    {return {true, true, false};}
                else
                    {return {false, false, false};}
            }
        }

    }

    // Section 2.3: Sufficient Conditions for eta. related and unrelated
    //// 2.3.1 Spacelike Bounds
    else if (r1 >= r2)
    {
        //spacelike
        if (t2-t1 < r1-r2)
            {return {false, false, false};}
        //spacelike
        else if (t2-t1 < r2*varphi2) 
            {return {false, false, false};}
        //timelike
        else if (r1 > 2*mass)
        {
            //First find r0
            double r0;
            if (r2 >= 3*mass){r0 = r2;}
            else if (r1 >= 3*mass && 3*mass > r2) {r0 = 3*mass;}
            else {r0 = r1;} //if (3*mass > r1 && r1 > 2*mass) 
            //then
            if (t2 >= t1 + r1 - r2 + (r0/std::sqrt(1-2*mass/r0))*varphi2)
                {return {true, true, false};}
            else
                {return BH_general_last_resort(transf_xvec, transf_yvec, mass);}
        }
        else
            {return BH_general_last_resort(transf_xvec, transf_yvec, mass);}
    }

    else if (r2 > r1 && r1 > 2*mass)
    {
        //spacelike
        if (t2-t1 < r2-r1 + 4*std::log((r2-2*mass)/(r1-2*mass)))
            {return {false, false, false};}
        else if (r1 > 3*mass)
        {
            //spacelike
            double r0 = r1;
            if (t2-t1 < (r0/std::sqrt(1-2*mass/r0))*varphi2)
                {return {false, false, false};}
            //timelike
            else if (t2>= t1 +r2 -r1 
                    + 4*mass * std::log((r2-2*mass)/(r1-2*mass))
                    + (r0/std::sqrt(1-2*mass/r0))*varphi2)
                {return {true, true, false};}
            else
                {return BH_general_last_resort(transf_xvec, transf_yvec, mass);}

        }
        else if(r1 < 3*mass && 3*mass < r2)
        {
            //spacelike
            double r0 = 3*mass;
            if (t2-t1 < (r0/std::sqrt(1-2*mass/r0))*varphi2)
                {return {false, false, false};}
            else
                {return BH_general_last_resort(transf_xvec, transf_yvec, mass);}
        }
        else if(r2 <= 3*mass)
        {
            //spacelike
            double r0 = r2;
            if (t2-t1 < (r0/std::sqrt(1-2*mass/r0))*varphi2)
                {return {false, false, false};}
            else
                {return BH_general_last_resort(transf_xvec, transf_yvec, mass);}
        }
        else
            {return BH_general_last_resort(transf_xvec, transf_yvec, mass);}
    }

    else //r2>2*mass>r1
        {return {false, false, false};}
}



/**
 * @brief Causality algorithm for two events in 4D EForig coordinates, from
 * Song He and David Rideout 2009 Class. Quantum Grav. 26 125015. 
 * 
 * @param xvec vector<double> : EF coordinates of x.
 * @param yvec vector<double> : EF coordinates of y.
 * @param period vector<double> : period along SPATIAL coordinates. Currently 
 * not implemented in BH, hence deafult is {}.
 * @param mass : mass of Black Hole
 * @return vector<bool> : {x-y timelike, x<=y, x>y}
 */
vector<bool> Spacetime::BH_general_causal4D (const vector<double>& xvec, 
                                            const vector<double>& yvec,
                                            std::vector<double> period,
                                            double mass)
{
    //WORKING IN EF COORDINATES
    if (yvec[0]<xvec[0])
    {
        std:vector<bool> result = Spacetime::BH_general_causal4D
                                             (yvec, xvec, period);
        if (result[0])
        {
            bool a = result[1]*1;
            result[1] = result[2]*1;
            result[2] = result[1]*1;
        }
        return result;
    }

    double t1     = xvec[0]; double t2     = yvec[0];
    double r1     = xvec[1]; double r2     = yvec[1];
    double theta1 = xvec[2]; double theta2 = yvec[2];
    double phi1   = xvec[3]; double phi2   = yvec[3];

    double vartheta1 = M_PI / 2;
    double vartheta2 = M_PI / 2; 
    double varphi1 = 0;
    double varphi2 = std::acos(std::cos(theta1)*std::cos(theta2) 
                    +std::sin(theta1)*std::sin(theta2)*std::cos(phi1-phi2));
    if (varphi2 < 0)
        {varphi2 += 2*M_PI;}
    
    vector<double> transf_xvec = {t1, r1, vartheta1, varphi1};
    vector<double> transf_yvec = {t2, r2, vartheta2, varphi2};
    
    // Section 2.2: Radially separated pairs and radial null geodesics
    if (varphi2<1e-6) //should be ==zero, but leave room for some error
    {
        if (r1>=r2)
        {
            // all 3 cases of the paper return same
            if (t2 >= t1 + r1 - r2 - 1e-6)
                    {return {true, true, false};}
            else
                {return {false, false, false};}
        }
        else /*r2>r1*/
        {
            if (r1<=2*mass)
                {return {false, false, false};}
            else // if (r1>2*mass)
            {
                if (t2 >= t1 + r2 - r1 + 
                            4*mass*std::log((r2-2*mass)/(r1-2*mass)))
                    {return {true, true, false};}
                else
                    {return {false, false, false};}
            }
        }

    }

    // Section 2.3: Sufficient Conditions for c. related and unrelated
    //// 2.3.1 Spacelike Bounds
    else if (r1 >= r2)
    {
        //spacelike
        if (t2-t1 < r1-r2)
            {return {false, false, false};}
        //spacelike
        else if (t2-t1 < r2*varphi2) 
            {return {false, false, false};}
        //timelike
        else if (r1 > 2*mass)
        {
            //First find r0
            double r0;
            if (r2 >= 3*mass){r0 = r2;}
            else if (r1 >= 3*mass && 3*mass > r2) {r0 = 3*mass;}
            else {r0 = r1;} //if (3*mass > r1 && r1 > 2*mass) 
            //then
            if (t2 >= t1 + r1 - r2 + (r0/std::sqrt(1-2*mass/r0))*varphi2)
                {return {true, true, false};}
            else
                {return BH_general_last_resort(transf_xvec, transf_yvec, mass);}
        }
        else
            {return BH_general_last_resort(transf_xvec, transf_yvec, mass);}
    }

    else if (r2 > r1 && r1 > 2*mass) /*r2>r1 and both outside horizon*/
    {
        //spacelike
        if (t2-t1 < r2-r1 + 4*std::log((r2-2*mass)/(r1-2*mass)))
            {return {false, false, false};}
        else if (r1 > 3*mass)
        {
            //spacelike
            double r0 = r1;
            if (t2-t1 < (r0/std::sqrt(1-2*mass/r0))*varphi2)
                {return {false, false, false};}
            //timelike
            else if (t2>= t1 +r2 -r1 
                    + 4*mass * std::log((r2-2*mass)/(r1-2*mass))
                    + (r0/std::sqrt(1-2*mass/r0))*varphi2)
                {return {true, true, false};}
            else
                {return BH_general_last_resort(transf_xvec, transf_yvec, mass);}

        }
        else if(r1 < 3*mass && 3*mass < r2)
        {
            //spacelike
            double r0 = 3*mass;
            if (t2-t1 < (r0/std::sqrt(1-2*mass/r0))*varphi2)
                {return {false, false, false};}
            else
                {return BH_general_last_resort(transf_xvec, transf_yvec, mass);}
        }
        else if(r2 <= 3*mass)
        {
            //spacelike
            double r0 = r2;
            if (t2-t1 < (r0/std::sqrt(1-2*mass/r0))*varphi2)
                {return {false, false, false};}
            else
                {return BH_general_last_resort(transf_xvec, transf_yvec, mass);}
        }
        else
            {return BH_general_last_resort(transf_xvec, transf_yvec, mass);}
    }

    else //r2>r1 with 2*mass>r1 
         //(i.e. r1 inside horizon, r2 could be inside or outside)
        {return {false, false, false};}
}


/**
 * @brief Last step in He, Rideout Algorithm for EF causality
 * 
 * @param xvec vector<double> : EF coordinates of x 
 * @param yvec vector<double> : EF coordinates of y
 * @param mass double : BH mass
 * @return vector<bool> : causality booleans
 */
vector<bool> Spacetime::BH_general_last_resort(const vector<double>& xvec, 
                                               const vector<double>& yvec,
                                               double mass)
{    
    double eta = Spacetime::BH_eta_solver(1./xvec[1],1./yvec[1],yvec[3], mass);
    if (eta<0)
    {
        return {false, false, false};
    }
    else
    {
        double geo_time = Spacetime::BH_int_dt_du (1./xvec[1],1./yvec[1], eta, 
                                                    mass);
        bool x_prec_y =  geo_time <= yvec[0] - xvec[0];
        return {x_prec_y, x_prec_y, false};
    }
}



// int main(){
// std::cout << "spacetimes.cpp WORKS! :)";
// }
