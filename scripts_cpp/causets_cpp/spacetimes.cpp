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



/*============================================================================
* ============================================================================
* GENERAL SPACETIME METHODS
* ============================================================================
============================================================================*/

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





/*============================================================================
* ============================================================================
* FLAT SPACETIME
* ============================================================================
============================================================================*/

/**
 * @brief Construct a new Flat Spacetime:: Flat Spacetime object
 * 
 * @param dim int:   dimension of flat spacetime. Default 4.
 * @param period vector<double>:   periodicity vector: ith entry is period 
 * along ith SPATIAL diemnsion (0->x, 1->y, 2->z). Default is no periodicty.
 * No periodicity along ith direction requires 0. Usually, but not necessarily,
 * with cuboid might want period[i]=shape_object.Edges()[i].
 */
void Spacetime::FlatSpacetime(int dim, vector<double> period)
{
    if (dim < 1)
    {
        std::cout<<"Given dim was smaller than 1."<<std::endl;
        throw std::invalid_argument("Dimension has to be at least 1.");
    }
    else if (period.size()!=0 && period.size()!=dim-1)
    {
        std::cout<<"Period size was "<<period.size()<<", it has to be either\
 0 or dim-1. Note: "<<"dimension was "<<dim<<std::endl;
        throw std::invalid_argument("Period's size is wrong!");
    }

    _dim = dim;
    _name = "Flat";
    _metricname = "Minkowski";
    if (period.size() && std::count(period.begin(), period.end(), 0)!=dim-1) 
    {
        _isPeriodic = true;
        _period = period;
        for (int i=0; i<period.size(); i++)
        {
            if (period[i]<0) {_period[i]=-period[i];}
        }
    }
}


/**
 * @brief Spacetime Interval Squared t squared - x_i*x^i
 */
double Spacetime::Flat_ds2(const vector<double>& xvec, const vector<double>& yvec)
{
    
    double time_delta2  = (yvec[0] - xvec[0])*(yvec[0] - xvec[0]);
    double space_delta2 = 0;
    for (int i = 1; i<_dim; i++)
        {space_delta2 += (yvec[i]-xvec[i])*(yvec[i]-xvec[i]);}
    return time_delta2 - space_delta2;
}


/**
 * @brief Spacetime Interval sqrt(t squared - x_i*x^i)
 */
double Spacetime::Flat_ds(const vector<double>& xvec, const vector<double>& yvec)
{
    
    double time_delta2  = (yvec[0] - xvec[0])*(yvec[0] - xvec[0]);
    double space_delta2 = 0;
    for (int i = 1; i<_dim; i++)
        {space_delta2 += (yvec[i]-xvec[i])*(yvec[i]-xvec[i]);}
    return std::sqrt(time_delta2 - space_delta2);
}


/**
 * @brief Function of two events in any D returning {x-y timelike?, x<=y, x>y}.
 * 
 * @param xvec vector<double>:  coordinates of x
 * @param yvec vector<double>:  coordinates of y
 * @param period not used, needed for consistency with Causality
 * @return bool : x-y timelike?
 */
bool Spacetime::Flat_causal(const vector<double>& xvec, const vector<double>& yvec,
                            vector<double>period, double mass)
{
    double dt = (xvec[0]-yvec[0]);
    double dspacex = xvec[1]-yvec[1];
    double dspacey = xvec[2]-yvec[2];
    double dspacez = xvec[3]-yvec[3];
    // for(int i=1; i<dim; i++){
    //    dspace2 += (xvec[i]-yvec[i])*(xvec[i]-yvec[i]);
    // }
    if ((dt*dt)-(dspacex*dspacex + dspacey*dspacey + dspacez*dspacez)>0)
    {
        return true;}
    else{
        return false;}
}

/**
 * @brief Function of two events in any D returning {x-y timelike?, x<=y, x>y}.
 * 
 * @param xvec : vector<double>, coordinates of x
 * @param yvec : vector<double>, coordinates of y
 * @param period : vector<double>, ith entry is periodicty along ith SPATIAL
 * dimension (0->x, 1->y, 2-->z)
 * @return bool : x-y timelike?
 */
bool Spacetime::Flat_causal_periodic(const vector<double>& xvec, 
                                        const vector<double>& yvec,
                                        vector<double> period,
                                        double mass)
{
    //std::cout<<"For debuggin: calling Flat_causal_periodic\n";
    double t_delta = (yvec[0] - xvec[0]);
    double t_delta2 = t_delta*t_delta;
    double space_delta2 = 0;
    for (int i = 1; i<yvec.size(); i++)
    {
        double space_delta_i = (yvec[i]-xvec[i]);
        double wrapped_space_delta_i = (space_delta_i>0)?
                                        period[i]-space_delta_i:
                                        period[i]+space_delta_i; 
        space_delta2 += (space_delta_i<wrapped_space_delta_i)?
                        space_delta_i*space_delta_i :
                        wrapped_space_delta_i*wrapped_space_delta_i;
    }
    return t_delta2 >= space_delta2;
}



/**
 * @brief Function of two events in any D returning {x-y timelike?, x<=y, x>y}.
 * 
 * @param xvec vector<double>:  coordinates of x
 * @param yvec vector<double>:  coordinates of y
 * @param period not used, needed for consistency with Causality
 * @return vector<bool> : {x-y timelike, x<=y, x>y}
 */
vector<bool> Spacetime::Flat_general_causal(const vector<double>& xvec, 
                                            const vector<double>& yvec,
                                            vector<double>period, double mass)
{
    double dt2 = (xvec[0]-yvec[0])*(xvec[0]-yvec[0]);
    double dspace2 = 0;
    for(int i=1; i<4; i++){
        dspace2 += (xvec[i]-yvec[i])*(xvec[i]-yvec[i]);
    }
    if (dt2-dspace2>0)
    {
        if (xvec[0]<yvec[0]) return {true,true, false};
        else return {true, false, true};
    }
    else
    {
        return {false, false, false};
    }
}


/**
 * @brief Function of two events in any D returning {x-y timelike?, x<=y, x>y}.
 * 
 * @param xvec : vector<double>, coordinates of x
 * @param yvec : vector<double>, coordinates of y
 * @param period : vector<double>, ith entry is periodicty along ith SPATIAL
 * dimension (0->x, 1->y, 2-->z)
 * @return vector<bool> : {x-y timelike, x<=y, x>y}
 */
vector<bool> Spacetime::Flat_general_causal_periodic(const vector<double>& xvec, 
                                                    const vector<double>& yvec,
                                                    vector<double> period,
                                                    double mass)
{
    //std::cout<<"For debuggin: calling Flat_causal_periodic\n";
    double t_delta = (yvec[0] - xvec[0]);
    double t_delta2 = t_delta*t_delta;
    double space_delta2 = 0;
    for (int i = 1; i<yvec.size(); i++)
    {
        double space_delta_i = (yvec[i]-xvec[i]);
        double wrapped_space_delta_i = (space_delta_i>0)?
                                        period[i]-space_delta_i:
                                        period[i]+space_delta_i; 
        space_delta2 += (space_delta_i<wrapped_space_delta_i)?
                        space_delta_i*space_delta_i :
                        wrapped_space_delta_i*wrapped_space_delta_i;
    }
    bool isCausal = t_delta2 >= space_delta2;
    return {isCausal,
            (t_delta >= 0.0) && isCausal,
            (t_delta < 0.0) && isCausal
           };
}







// /*============================================================================
// * ============================================================================
// * BLACK HOLES SPACETIME
// * Notes: does not support anything actually at the moment
// * ============================================================================
// ============================================================================*/

/**
 * @brief Initialises a Black Hole Spacetime.
 *
 * @param dim: dimension of spacetime. Default 2.
 * @param mass: Schwarzschild mass. Default 1.
 * @param metric: Specify metric: either "EF(original)" (default), or "EF(uv)",
 *                or "Schwarzschild" or "S".
 */
void Spacetime::BlackHoleSpacetime(int dim,// = 4
                                    double mass,// = 1
                                    std::string metric)// = "EF(original)"
{

    if (dim != 4 && dim != 2 && dim !=3)
    {
        std::cout<<"Dimension has to be 2, 3 or 4."<<std::endl;
        throw std::invalid_argument("Dimension has to be 2,3 or 4.");
    }
    _dim = dim;
    _name = "BlackHole";
    _mass = mass;
    _r_S  = 2*mass;

    if (metric == "Eddington-Finkelstein(original)" || metric=="EF(original)")
        {_metricname = "EF(original)";}
    else if (metric == "Eddington-Finkelstein(uv)" || metric == "EF(uv)")
        {_metricname = "EF(uv)";}
    else if (metric == "Schwarzschild" || metric == "S")
        {_metricname = "Schwarzschild";}
    else 
    {
        std::cout<<"Given metric for BlackHole must be 'Schwarzschild', "<<
                 "'EF(uv)' or 'EF(original)'. "<<std::endl;
        throw std::invalid_argument("Wrong metricname");
    }
}


/**
 * @brief Causality algorithm for two events in 2D EForig coordinates, from
 * Song He and David Rideout 2009 Class. Quantum Grav. 26 125015. 
 * 
 * @param xvec vector<double> : EF coordinates of x.
 * @param yvec vector<double> : EF coordinates of y.
 * @param period vector<double> : period along SPATIAL coordinates. Currently 
 * not implemented in BH, hence deafult is {}.
 * @param mass : mass of Black Hole
 * @return bool : x-y timelike?
 */
bool Spacetime::BH_causal2D (const vector<double>& xvec, 
                             const vector<double>& yvec,
                             std::vector<double> period,
                             double mass)
{
    //IF WORKING IN EF COORDINATES

    if (yvec[0]<xvec[0])
        {return Spacetime::BH_causal3D(yvec, xvec, period);}

    double t1     = xvec[0]; double t2     = yvec[0];
    double r1     = xvec[1]; double r2     = yvec[1];

    if (r1*r2 < 0) //on opposite sides of the singularity
        {return false;}
    else
    {
        r1 = std::abs(r1);
        r2 = std::abs(r2);
    }
    
    // Section 2.2: Radially separated pairs and radial null geodesics
    // as in 2D necessarily varphi2 = 0
    if (r1>=r2)
    {
        // all 3 cases of the paper require this condition
        if (t2 >= t1 + r1 - r2)
        {
            //the further is outside
            if (r1> 2*mass)
                {return true;}
           
            //both inside
            else
                {
                    if (t2 <= t1 + r2 - r1
                              + 4*mass*std::log( (2*mass-r2)/(2*mass-r1) )
                        )
                        return true; 
                    else 
                        return false;
                }
        }
        else
            {return false;}
    }
    else /*r2>r1*/
    {
        if (r1<=2*mass)
            {return false;}
        else // if (r1>2*mass)
        {
            if (t2 >= t1 + r2 - r1 + 4*mass*std::log((r2-2*mass)/(r1-2*mass)))
                {return true;}
            else
                {return false;}
        }
    }
}


/**
 * @brief Causality algorithm for two events in 4D EForig coordinates, mainly 
 * from Song He and David Rideout 2009 Class. Quantum Grav. 26 125015. However,
 * in their paper it was assumed any object could move on a constant-spatial-
 * components line everywhere. This doesn't hold inside the horizon, therefore,
 * a new upper bound was required.
 * 
 * @param xvec vector<double> : EF coordinates of x.
 * @param yvec vector<double> : EF coordinates of y.
 * @param period vector<double> : period along SPATIAL coordinates. Currently 
 * not implemented in BH, hence deafult is {}.
 * @param mass : mass of Black Hole
 * 
 * @return bool : x-y timelike?
 */
bool Spacetime::BH_causal3D (const std::vector<double>& xvec, 
                            const std::vector<double>& yvec,
                            std::vector<double> period,
                            double mass)
{
    if (yvec[0]<xvec[0])
        {return Spacetime::BH_causal3D(yvec, xvec, period);}
    
    else
    {
        const std::vector<double> xvec4D = {xvec[0], xvec[1], M_PI/2., xvec[2]};
        const std::vector<double> yvec4D = {yvec[0], yvec[1], M_PI/2., yvec[2]};
        return BH_causal4D(xvec4D, yvec4D, period, mass);
    }

    // double t1     = xvec[0]; double t2     = yvec[0];
    // double r1     = xvec[1]; double r2     = yvec[1];
    // double phi1   = xvec[2]; double phi2   = yvec[2];
    //
    // //double vartheta1 = M_PI / 2;
    // //double vartheta2 = M_PI / 2; 
    // //double varphi1 = 0;
    // double varphi2 = phi1-phi2;
    // if (varphi2 < 0)
    //     {varphi2 += 2*M_PI;}
    //
    // vector<double> transf_xvec = {t1, r1, M_PI / 2, 0}; //{t1, r1, vartheta1, varphi1};
    // vector<double> transf_yvec = {t2, r2, M_PI / 2, varphi2}; //{t2, r2, vartheta2, varphi2};
    //
    // // Section 2.2: Radially separated pairs and radial null geodesics
    // if (varphi2<1e-6) //should be ==zero, but leave room for some error
    // {
    //     if (r1>=r2)
    //     {
    //         // all 3 cases of the paper require this condition
    //         if (t2 >= t1 + r1 - r2)
    //         {
    //             //the further is outside
    //             if (r1> 2*mass)
    //                 {return true;}
    //        
    //             //both inside
    //             else
    //                 {
    //                     if (t2 <= t1 + r2 - r1
    //                             + 4*mass*std::log( (2*mass-r2)/(2*mass-r1) )
    //                         )
    //                         return true; 
    //                     else 
    //                         return false;
    //                 }
    //         }
    //         else
    //             {return false;}
    //     }
    // }
    //
    // // Section 2.3: Sufficient Conditions for c. related and unrelated
    // //// 2.3.1 Spacelike Bounds
    // else if (r1 >= r2)
    // {
    //     //spacelike
    //     if (t2-t1 < r1-r2)
    //         {return false;}
    //     //spacelike
    //     else if (t2-t1 < r2*varphi2) 
    //         {return false;}
    //     //timelike
    //     else if (r1 > 2*mass)
    //     {
    //         //First find r0
    //         double r0;
    //         if (r2 >= 3*mass){r0 = r2;}
    //         else if (r1 >= 3*mass && 3*mass > r2) {r0 = 3*mass;}
    //         else {r0 = r1;} //if (3*mass > r1 && r1 > 2*mass) 
    //         //then
    //         if (t2 >= t1 + r1 - r2 + (r0/std::sqrt(1-2*mass/r0))*varphi2)
    //             {return true;}
    //         else
    //             {return BH_last_resort(transf_xvec, transf_yvec, mass);}
    //     }
    //     else
    //         {return BH_last_resort(transf_xvec, transf_yvec, mass);}
    // }
    //
    // else if (r2 > r1 && r1 > 2*mass)
    // {
    //     //spacelike
    //     if (t2-t1 < r2-r1 + 4*std::log((r2-2*mass)/(r1-2*mass)))
    //         {return false;}
    //     else if (r1 > 3*mass)
    //     {
    //         //spacelike
    //         double r0 = r1;
    //         if (t2-t1 < (r0/std::sqrt(1-2*mass/r0))*varphi2)
    //             {return false;}
    //         //timelike
    //         else if (t2>= t1 +r2 -r1 
    //                 + 4*mass * std::log((r2-2*mass)/(r1-2*mass))
    //                 + (r0/std::sqrt(1-2*mass/r0))*varphi2)
    //             {return true;}
    //         else
    //             {return BH_last_resort(transf_xvec, transf_yvec, mass);}
    //
    //     }
    //     else if(r1 < 3*mass && 3*mass < r2)
    //     {
    //         //spacelike
    //         double r0 = 3*mass;
    //         if (t2-t1 < (r0/std::sqrt(1-2*mass/r0))*varphi2)
    //             {return false;}
    //         else
    //             {return BH_last_resort(transf_xvec, transf_yvec, mass);}
    //     }
    //     else if(r2 <= 3*mass)
    //     {
    //         //spacelike
    //         double r0 = r2;
    //         if (t2-t1 < (r0/std::sqrt(1-2*mass/r0))*varphi2)
    //             {return false;}
    //         else
    //             {return BH_last_resort(transf_xvec, transf_yvec, mass);}
    //     }
    //     else
    //         {return BH_last_resort(transf_xvec, transf_yvec, mass);}
    // }
    //
    // else //r2>2*mass>r1
    //     {return false;}
    //
    // return false;
}



/**
 * @brief Causality algorithm for two events in 4D EForig coordinates, mainly 
 * from Song He and David Rideout 2009 Class. Quantum Grav. 26 125015. However,
 * in their paper it was assumed any object could move on a constant-spatial-
 * components line everywhere. This doesn't hold inside the horizon, therefore,
 * a new upper bound was required. 
 * 
 * @param xvec vector<double> : EF coordinates of x.
 * @param yvec vector<double> : EF coordinates of y.
 * @param period vector<double> : period along SPATIAL coordinates. Currently 
 * not implemented in BH, hence deafult is {}.
 * @param mass : mass of Black Hole
 * @return bool : x-y timelike?
 */
bool Spacetime::BH_causal4D (const vector<double>& xvec, 
                             const vector<double>& yvec,
                             std::vector<double> period,
                             double mass)
{
    double t2_min_t1 = yvec[0] - xvec[0];
    double r2_min_r1 = yvec[1] - xvec[1];
    double varphi2 = std::acos(std::cos(xvec[2])*std::cos(yvec[2]) 
            +std::sin(xvec[2])*std::sin(yvec[2])*std::cos(xvec[3]-yvec[3]));
    if (varphi2 < 0)
        {varphi2 += 2*M_PI;}
    
    vector<double> transf_xvec = {xvec[0], xvec[1], M_PI/2., 0};
    vector<double> transf_yvec = {yvec[0], yvec[1], M_PI/2., varphi2};
    
    // Section 2.2: Radially separated pairs and radial null geodesics
    if (varphi2<1e-6) //should be ==zero, but leave room for some error
    {
        if (r2_min_r1<=0) /*r1>=r2*/
        {
            // all 3 cases of the paper require this condition
            if (t2_min_t1 + r2_min_r1 >= 0) /*t2 >= t1 + r1 - r2*/
            {
                //the further is outside
                if (xvec[1] > 2*mass) /*r1 > 2M*/
                    {return true;}
            
                //both inside
                else
                    {
                        if (t2_min_t1 <= r2_min_r1
                                + 4*mass*std::log( (2*mass-yvec[1])
                                                 / (2*mass-xvec[1]) ))
                            return true; 
                        else 
                            return false;
                    }
            }
            else
                {return false;}
        }
        else
        {
            if (xvec[1] - 2*mass<=0)
                {return false;}
            else // if (r1>2*mass)
            {
                if (t2_min_t1 - r2_min_r1 - 
                    4*mass*std::log((yvec[1]-2*mass)/(xvec[1]-2*mass)) >= 0)
                    {return true;}
                else
                    {return false;}
            }
        }
    }

    // Section 2.3: Sufficient Conditions for c. related and unrelated
    //// 2.3.1 Spacelike Bounds
    else if (r2_min_r1 <= 0)// xvec[1] >= yvec[1])
    {
        //spacelike
        if (t2_min_t1 + r2_min_r1 < 0)
            {return false;}
        //spacelike
        else if (t2_min_t1 - yvec[1]*varphi2 < 0) 
            {return false;}
        //timelike
        else if (xvec[1] - 2*mass> 0)
        {
            //First find r0
            double r0;
            if (yvec[1] - 3*mass>= 0){r0 = yvec[1];}
            else if (xvec[1] >= 3*mass && 3*mass > yvec[1]) {r0 = 3*mass;}
            else {r0 = xvec[1];} //if (3*mass > r1 && r1 > 2*mass) 
            //then
            if (t2_min_t1 + r2_min_r1 >= (r0/std::sqrt(1-2*mass/r0))*varphi2)
                {return true;}
            else
                {return BH_last_resort(transf_xvec, transf_yvec, mass);}
        }
        else
            {return BH_last_resort(transf_xvec, transf_yvec, mass);}
    }

    else if (r2_min_r1>0 && xvec[1] - 2*mass> 0)
    {
        //spacelike
        if (t2_min_t1-r2_min_r1-4*std::log((yvec[1]-2*mass)/(xvec[1]-2*mass)) < 0)
            {return false;}
        else if (xvec[1] - 3*mass> 0)
        {
            //spacelike
            double r0 = xvec[1];
            if (t2_min_t1 < (r0/std::sqrt(1-2*mass/r0))*varphi2)
                {return false;}
            //timelike
            else if (t2_min_t1>= r2_min_r1 
                    + 4*mass * std::log((yvec[1]-2*mass)/(xvec[1]-2*mass))
                    + (r0/std::sqrt(1-2*mass/r0))*varphi2)
                {return true;}
            else
                {return BH_last_resort(transf_xvec, transf_yvec, mass);}

        }
        else if(xvec[1] - 3*mass< 0 && 3*mass - yvec[1] < 0)
        {
            //spacelike
            double r0 = 3*mass;
            if (t2_min_t1 < (r0/std::sqrt(1-2*mass/r0))*varphi2)
                {return false;}
            else
                {return BH_last_resort(transf_xvec, transf_yvec, mass);}
        }
        else if(yvec[1] - 3*mass<= 0)
        {
            //spacelike
            double r0 = yvec[1];
            if (t2_min_t1 < (r0/std::sqrt(1-2*mass/r0))*varphi2)
                {return false;}
            else
                {return BH_last_resort(transf_xvec, transf_yvec, mass);}
        }
        else
            {return BH_last_resort(transf_xvec, transf_yvec, mass);}
    }

    else //r2>2*mass>r1
        {return false;}
}


/**
 * @brief Integral last step in He, Rideout Algorithm for EForig causality. 
 * An upper bound was added to account for the impossibility of a stationary
 * object inside the Horizon. Uses eta = E_S/L, not eta = E^{*}/L 
 * (which would have its sign inverted).
 * 
 * @param xvec vector<double> : EF coordinates of x 
 * @param yvec vector<double> : EF coordinates of y
 * @param mass double : BH mass
 * 
 * @return bool : causality booleans
 */
bool Spacetime::BH_last_resort(const vector<double>& xvec, 
                                const vector<double>& yvec,
                                double mass)
{
    double eta = Spacetime::BH_eta_solver(1./xvec[1],1./yvec[1],yvec[3], mass);
    //no sol found
    if (eta<0)
    {
        //std::cout<<"   eta^2          is :"<<-1      <<std::endl;
        return false;
    }
    else
    {
        // note as eta = E_S/L rather than E^{*}/L, positive eta gives lower
        // bound
        double geo_time1 = Spacetime::BH_int_dt_du (1./xvec[1],1./yvec[1], eta, 
                                                    mass);
        bool x_prec_y;
        if (xvec[1] < 2*mass && yvec[1] < 2*mass) /*both inside*/
        {
            double geo_time2 = Spacetime::BH_int_dt_du (1./xvec[1],1./yvec[1],
                                                         -eta, mass);
            x_prec_y =  geo_time1 <= yvec[0] - xvec[0]
                       && yvec[0] - xvec[0] <= geo_time2;
        }
        else
        {
            x_prec_y =  geo_time1 <= yvec[0] - xvec[0];
        }
        //std::cout<<"   eta^2          is :"<<eta*eta      <<std::endl;
        //std::cout<<"geodesic's time is : "<<geo_time<<std::endl;
        return x_prec_y;
    }
}


/**
 * @brief d(varphi)/du for a BH as of He and Rideout (Eq.13).
 * 
 * @param dpdu double& : derivative gets updated.
 * @param u double : u = 1/r.
 * @param eta2 double : eta constant squared.
 * @param M double : mass of BH.
 */
void Spacetime::BH_dvarphi_du (double& dpdu, double u, double eta2, double M)
    {dpdu = std::pow(2*M*u*u*u - u*u + eta2, -0.5);}


/**
 * @brief Integral of d(varphi)/du from u1 to u2.
 * 
 * @param u1 double : 1/r_1
 * @param u2 double : 1/r_2
 * @param eta2 double : eta constant squared.
 * @param M double : mass of BH.
 * @return double : the result of the integral. 
 */
double Spacetime::BH_int_dvarphi_du(double u1, double u2, double eta2, double M)
{
    auto BH_dvarphi_du_forint = [M, eta2]
                            (const double& varphi, double& dpdu, const double u)
                            {dpdu = std::pow(u*u*(2*M*u - 1) + eta2, -0.5);};
    double varphi = 0;
    if (u2>=u1)
    {
        boost::numeric::odeint::integrate(BH_dvarphi_du_forint, varphi, 
                                          u1, u2, (u2-u1)/20.);
    }
    else
    {
        boost::numeric::odeint::integrate(BH_dvarphi_du_forint, varphi, 
                                          u2, u1, (u1-u2)/20.);
    }
    return varphi;
}


/**
 * @brief Fit the integral of d(varphi)/du to varphi2 to get best estimate of 
 * eta^2 by bisection method (where eta=E/L)
 * 
 * @param u1 double : 1/r_1
 * @param u2 double : 1/r_2
 * @param varphi2 double : the value we want the integral to be.
 * @param M double : mass of BH.
 * 
 * @return double : the best fit for eta=E/L, up to uncertainty 1e-3.
 */
double Spacetime::BH_eta_solver (double u1, double u2, double varphi2, double M)
{
  // U2>=U1 -> USE PLUS
  if (u2>=u1) 
  {
    //SET MINIMUM BOUND
    double eta2min = 0;
    if (2*M*u1<=1.)
    {
        double D1 = u1*u1 *(1 - 2*M*u1);
        double D2 = u2*u2 *(1 - 2*M*u2);
        eta2min += (D1>D2)? D1 : D2;
        eta2min += 1e-5; //for anti-divergence purposes
    }
    //SET UPPER BOUND
    double deltaphi_max = BH_int_dvarphi_du(u1, u2, eta2min, M);
    if (deltaphi_max <= varphi2)
        {return -1;}//std::sqrt(eta2min);}
    else
    {
        double eta = (u2-u1)/varphi2;
        double eta2max = 0;
        if (2*M*u1 >= 1.)
            {eta2max += eta*eta;}
        else
            {eta2max += (u2+eta)*(u2+eta);}
        // check upper bound is on other side of solution. Since deltaphi_max
        // is positive, this has to be negative. If not, double upper limit,
        // turn lower limit to previous upper and check again.
        while (BH_int_dvarphi_du(u1, u2, eta2max, M) > varphi2)
        {
            std::cout<<"No biggy, just in c_solver had to update upper limit"
                    <<std::endl;
            eta2min = eta2max*1;
            eta2max *= 2;
        }
        // SOLVE WITH BISECTION
        auto BH_tosolve = [u1, u2, varphi2, M](double eta2)
                        {return BH_int_dvarphi_du(u1, u2, eta2, M)-varphi2;};
        return std::sqrt(bisection(BH_tosolve, eta2min, eta2max, 1e-5));
    }
  }
  // U2<U1 -> USE MINUS (in dvarphi_du it means switch u1 and u2)
  else 
  {
        //SET MINIMUM BOUND
        double eta2min = 0;
        if (2*M*u2<=1.)
        {
            double D1 = u1*u1 *(1 - 2*M*u1);
            double D2 = u2*u2 *(1 - 2*M*u2);
            eta2min += (D1>D2)? D1 : D2;
            eta2min += 1e-5; //for anti-divergence purposes
        }
        //SET UPPER BOUND
        double deltaphi_max = BH_int_dvarphi_du(u2, u1, eta2min, M);
        if (deltaphi_max <= 0)
        {return -1;}//std::sqrt(eta2min);}
        else
        {
        double eta = (u1-u2)/varphi2;
        double eta2max = 0;
        if (2*M*u2 >= 1)
            {eta2max += eta*eta;}
        else
            {eta2max += (u1+eta)*(u1+eta);}
        // check upper bound is on other side of solution. Since deltaphi_max
        // is positive, this has to be negative. If not, switch limits and check
        // again.
        while (BH_int_dvarphi_du(u2, u1, eta2max, M) > varphi2)
        {
            //std::cout<<"No biggy, just in c_solver had to update upper limit"
                //       <<std::endl;
            eta2min = eta2max*1;
            eta2max *= 2;
        }
        // SOLVE WITH BISECTION
        auto BH_tosolve = [u2, u1, varphi2, M](double eta2)
                        {return BH_int_dvarphi_du(u2, u1, eta2, M)-varphi2;};
        return std::sqrt(bisection(BH_tosolve, eta2min, eta2max, 1e-5));
    }
  }
}


/**
 * @brief dt/du for a BH as of He and Rideout (Eq.15), 
 * for u2>u1 in the integral.
 * 
 * @param dpdu double& : derivative gets updated.
 * @param u double : u = 1/r.
 * @param eta double : eta constant E/L (NOT squared)
 * @param M double : mass of BH.
 */
void Spacetime::BH_dt_du_plus (double&dtdu, double u, double eta, double M)
{
  double D = u*u*(1-2*M*u);
  dtdu = (eta*std::pow(eta*eta - D, -0.5) - 2*M*u)/D;
}


/**
 * @brief dt/du for a BH as of He and Rideout (Eq.15), 
 * for u1>u2 in the integral.
 * 
 * @param dpdu double& : derivative gets updated.
 * @param u double : u = 1/r.
 * @param eta double : eta constant E/L (NOT squared)
 * @param M double : mass of BH.
 */
void Spacetime::BH_dt_du_minus (double&dtdu, double u, double eta, double M)
{
  double D = u*u*(1-2*M*u);
  dtdu = (-eta*std::pow(eta*eta - D, -0.5) - 2*M*u)/D;
}


/**
 * @brief Integral of dt/du for a BH as of He and Rideout (Eq.15). Uses 
 * eta = E_S/L, not eta = E^{*}/L (which would have its sign inverted)
 * 
 * @param u1 double : u=1/r1
 * @param u2 double : u = 1/r2.
 * @param eta double : eta constant E/L (NOT squared)
 * @param M double : mass of BH.
 */
double Spacetime::BH_int_dt_du (double u1, double u2, double eta, double M)
{
    double t = 0;
    if (u2>=u1) //du>0
    {
        auto BH_dt_du_forint_plus = [M, eta]
                                (const double& t, double& dtdu, const double u)
                                {
                                  if (u==0.5)
                                    {dtdu = 0.488889;}
                                  else
                                  {
                                    double D = u*u*(1-2*M*u);
                                    dtdu = (eta/std::sqrt(eta*eta - D) - 2*M*u)
                                          /D;
                                  }
                                };
        //std::cout << "u2>=u1\n";
        boost::numeric::odeint::integrate(BH_dt_du_forint_plus, t, 
                                              u1, u2, (u2-u1)/20.);
    }
    else /* u1>u2, //du<0 */
    {
        auto BH_dt_du_forint_minus = [M, eta]
                                (const double& t, double& dtdu, const double u)
                                {
                                  double D = u*u*(1-2*M*u);
                                  dtdu = (-eta/std::sqrt(eta*eta - D) - 2*M*u)
                                        /D;
                                };
        if (0.5*M>=u1 || u2>=0.5*M) //0.5*M NOT in [u2, u1]
        {
            //std::cout << "in 0.5*M NOT in [u2, u1]\n";
            //Avoid divergence at 1-2*M*u=0
            u2 += (u2==0.5*M)? 1e-3*M : 0;
            u1 -= (u2==0.5*M)? 1e-3*M : 0;
            //Compute 
            boost::numeric::odeint::integrate(BH_dt_du_forint_minus, t, 
                                              u1, u2, -(u1-u2)/20.);
        }
        else /*0.5*M IN [u2, u1]*/
        {
            //std::cout << "in else.................................... \n";
            //Compute in 2 steps to avoid divergence
            
            boost::numeric::odeint::integrate(BH_dt_du_forint_minus, t, 
                                            u1, (0.5+0.0001)*M, -(u1-0.5)/20.);
            boost::numeric::odeint::integrate(BH_dt_du_forint_minus, t, 
                                            (0.5-0.0001)*M, u2, -(0.5-u2)/20.);
        }
    }
    return t;
}


bool Spacetime::BH_time_caus_check(double u1, double u2, double t1, double t2,
                                     double eta, double M)
    {return Spacetime::BH_int_dt_du (u1, u2, eta, M) + t1 - t2 <= 0;}




///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
// BH Coordinate Transformations
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////


typedef void (*inversefunc)
(std::vector<std::vector<double>>& coords, double mass, const char* EFtype);
/**
 * @brief Turn coords from "Spacetime_object._metricname" to "EF(original)" 
 * 
 * @param coords vector<vector<double>>& : coordinates to change
 * @return func "void callable(vector<vector<double>> &coords)" : function that
 * turns coordinates back from "EF(original)" to the first type.
 */
inversefunc Spacetime::ToInEF_original(std::vector<std::vector<double>>&coords)
{
    if (_metricname=="EF(uv)")
    {
        Spacetime::switchInEF(coords, "uv");
        return Spacetime::EF_from_original_to_uv;
    }
    else if (_metricname=="Schwarzschild")
    {
        Spacetime::StoInEF(coords, _mass, "original");
        return Spacetime::InEFtoS;
    }
    else if (_metricname=="GP")
    {
        Spacetime::GPtoInEF(coords, _mass);
        return Spacetime::InEFtoGP;
    }
    else
        {return Spacetime::do_nothing;}
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
        xvec[2] = phi;
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
