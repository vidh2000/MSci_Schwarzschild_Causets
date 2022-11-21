/// \authors Vid Homsak, Stefano Veroni
/// \date 25/09/2022

#ifndef SPACETIMES_H
#define SPACETIMES_H

#include <algorithm>
#include <cmath>
#include <cstring>
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

#include "boost/numeric/odeint.hpp"
//#include "shapes.h"


class Spacetime
{
    public:
        int  _dim;
        const char* _name;      //"Flat", "BlackHole"
        std::string _metricname;//"Minkowski", "Eddington-Finkelstein", etc... 

        //For BH
        double _mass = 0; 
        double _r_S = 0;

        //For Periodicity
        bool _isPeriodic = false;
        std::vector<double> _period = {};


        // CONSTRUCTOR

        Spacetime();


        // GENERAL FUNCTIONS

        std::vector<double>T_slice_sampling(double t, 
                                            std::vector<double>origin,
                                            int samplingsize = -1); 


        // CAUSALITY
        typedef std::vector<bool> (*func)
        (std::vector<double> xvec, std::vector<double> yvec, 
         std::vector<double> period, double mass);
        func Causality();  
        
        static std::vector<bool> causal1d(std::vector<double> xvec, 
                                          std::vector<double> yvec,
                                          std::vector<double> period,
                                          double mass);


    /*==============================================================
    ===============================================================*/
    //class FlatSpacetime: public Spacetime
//{
    public:
        
        // Constructor-like method

        void FlatSpacetime(int dim = 4, std::vector<double> period = {});

        // Flat Spacetime methods

        double Flat_ds2    (std::vector<double> xvec, std::vector<double> yvec);
        double Flat_ds     (std::vector<double> xvec, std::vector<double> yvec);  

        static std::vector<bool> Flat_causal (std::vector<double> xvec, 
                                         std::vector<double> yvec,
                                         std::vector<double> period = {},
                                         double mass = 0);
        static std::vector<bool> Flat_causal_periodic (std::vector<double> xvec, 
                                            std::vector<double> yvec,
                                            std::vector<double>period,
                                            double mass = 0);


    /*==============================================================
    ===============================================================*/
    //class BlackHoleSpacetime: public Spacetime

    public:

        // Constructor-like method

        void BlackHoleSpacetime(int dim = 4,
                                double mass = 1,
                                std::string metric = "Eddington-Finkelstein");
    
        // BH Causality

        double BH_ds2(std::vector<double> xvec, std::vector<double> yvec);
        double BH_ds(std::vector<double> xvec, std::vector<double> yvec); 

        static std::vector<bool> BH_causal2D (std::vector<double> xvec, 
                                              std::vector<double> yvec,
                                              std::vector<double> period={},
                                              double mass = 0);

        static std::vector<bool> BH_causal3D (std::vector<double> xvec, 
                                              std::vector<double> yvec,
                                              std::vector<double> period={},
                                              double mass = 0); 

        static std::vector<bool> BH_causal4D (std::vector<double> xvec, 
                                              std::vector<double> yvec,
                                              std::vector<double> period={},
                                              double mass = 0);
        
        static std::vector<bool> BH_last_resort(std::vector<double> xvec, 
                                                std::vector<double> yvec,
                                                double mass = 0);
        
        static void BH_dvarphi_du (double& dpdu, double u, double c2, 
                                   double M);
        static double BH_int_dvarphi_du(double u1, double u2, double c2,
                                        double M);
        static double BH_c_solver (double u1, double u2, double varphi2, 
                                   double M);

        static void BH_dt_du_plus  (double&dtdu, double u, double c,
                                    double M);
        static void BH_dt_du_minus (double&dtdu, double u, double c,
                                    double M);
        static double BH_int_dt_du (double u1, double u2, double c,
                                    double M);
        static bool BH_time_caus_check(double u1,double u2,double t1,double t2,
                                       double c, double M);


        // BH Coordinate Transformations

        static void InEFtoS (std::vector<double>& xvec, double mass = 1,
                             const char* EFtype = "original");
        static void InEFtoS (std::vector<std::vector<double>>& coords, 
                             double mass = 1,
                             const char* EFtype = "original");
        static void StoInEF (std::vector<double>& xvec, double mass = 1,
                             const char* EFtype = "original");
        static void StoInEF (std::vector<std::vector<double>>& coords,
                             double mass = 1,
                             const char* EFtype = "original");

        static void GPtoS (std::vector<double>& xvec, double mass = 1);
        static void GPtoS (std::vector<std::vector<double>>& coords, 
                             double mass = 1);
        static void StoGP (std::vector<double>& xvec, double mass = 1);
        static void StoGP (std::vector<std::vector<double>>& coords, 
                             double mass = 1);
        
        static void InEFtoGP (std::vector<double>& xvec, double mass = 1,
                             const char* EFtype = "original");
        static void InEFtoGP (std::vector<std::vector<double>>& coords, 
                             double mass = 1,
                             const char* EFtype = "original");
        static void GPtoInEF (std::vector<double>& xvec, double mass = 1,
                             const char* EFtype = "original");
        static void GPtoInEF (std::vector<std::vector<double>>& coords,
                             double mass = 1,
                             const char* EFtype = "original");
};

#endif /* SPACETIMES_H */