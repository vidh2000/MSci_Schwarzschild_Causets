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
        std::string _metricname;//"Minkowski",
                                // "Schwarzschild","EF(original),"EF(uv)","GP" 

        //For BH, but need them in general
        double _mass = 0; 
        double _r_S = 0;

        //For Periodicity
        bool _isPeriodic = false;
        std::vector<double> _period = {};


        Spacetime();

        std::vector<double>T_slice_sampling(double t, 
                                            std::vector<double>origin,
                                            int samplingsize = -1); 


        // CAUSALITY
        typedef bool (*func)
        (const std::vector<double>& xvec, const std::vector<double>& yvec, 
         std::vector<double> period, double mass);
        func Causality();

        typedef std::vector<bool> (*general_func)
        (const std::vector<double>& xvec, const std::vector<double>& yvec, 
         std::vector<double> period, double mass);
        general_func General_Causality();  
        
        static bool              causal1d(const std::vector<double>& xvec, 
                                          const std::vector<double>& yvec,
                                          std::vector<double> period,
                                          double mass);
        
        static std::vector<bool> general_causal1d(const std::vector<double>& xvec, 
                                                  const std::vector<double>& yvec,
                                                  std::vector<double> period,
                                                  double mass);
        

        static void CarttoSpherical (std::vector<double>& xvec);
        static void CarttoSpherical (std::vector<std::vector<double>>& coords);
        static void SphericaltoCart (std::vector<std::vector<double>>& coords);


    /*==============================================================
    ===============================================================*/
    //class FlatSpacetime: public Spacetime
//{
    public:
        
        // Constructor-like method

        void FlatSpacetime(int dim = 4, std::vector<double> period = {});

        // Flat Spacetime methods

        double Flat_ds2    (const std::vector<double>& xvec, 
                            const std::vector<double>& yvec);
        double Flat_ds     (const std::vector<double>& xvec, 
                            const std::vector<double>& yvec);  

        static bool Flat_causal (const std::vector<double>& xvec, 
                                 const std::vector<double>& yvec,
                                 std::vector<double> period = {},
                                 double mass = 0);
        static bool Flat_causal_periodic (const std::vector<double>& xvec, 
                                          const std::vector<double>& yvec,
                                          std::vector<double>period,
                                          double mass = 0);
        
        static std::vector<bool> Flat_general_causal (
                                            const std::vector<double>& xvec, 
                                            const std::vector<double>& yvec,
                                            std::vector<double> period = {},
                                            double mass = 0);
        static std::vector<bool> Flat_general_causal_periodic (
                                            const std::vector<double>& xvec, 
                                            const std::vector<double>& yvec,
                                            std::vector<double>period,
                                            double mass = 0);


    /*==============================================================
    ===============================================================*/
    //class BlackHoleSpacetime: public Spacetime

    public:

        // Constructor-like method

        void BlackHoleSpacetime(int dim = 4,
                                double mass = 1,
                                std::string metric = "EF(original)");
    
        // BH Causality

        double BH_ds2(const std::vector<double>& xvec, const std::vector<double>& yvec);
        double BH_ds(const std::vector<double>& xvec, const std::vector<double>& yvec); 

        static bool BH_causal2D (const std::vector<double>& xvec, 
                                 const std::vector<double>& yvec,
                                 std::vector<double> period={},
                                 double mass = 1);

        static bool BH_causal3D (const std::vector<double>& xvec, 
                                 const std::vector<double>& yvec,
                                 std::vector<double> period={},
                                 double mass = 1); 

        static bool BH_causal4D (const std::vector<double>& xvec, 
                                 const std::vector<double>& yvec,
                                 std::vector<double> period={},
                                 double mass = 1);
        
        static bool BH_last_resort(const std::vector<double>& xvec, 
                                   const std::vector<double>& yvec,
                                   double mass = 1);
        
        static void BH_dvarphi_du (double& dpdu, double u, double c2, 
                                   double M);
        static double BH_int_dvarphi_du(double u1, double u2, double c2,
                                        double M);
        static double BH_eta_solver (double u1, double u2, double varphi2, 
                                   double M);

        static void BH_dt_du_plus  (double&dtdu, double u, double c,
                                    double M);
        static void BH_dt_du_minus (double&dtdu, double u, double c,
                                    double M);
        static double BH_int_dt_du (double u1, double u2, double c,
                                    double M);
        static bool BH_time_caus_check(double u1,double u2,double t1,double t2,
                                       double c, double M);


        // BH Coordinate Transformations (GP NOT TESTED)

        typedef void (*inversefunc)
        (std::vector<std::vector<double>>& coords, double mass, 
         const char* EFtype);
        inversefunc ToInEF_original(std::vector<std::vector<double>> &coords);

        static void CarttoSpherical (std::vector<double>& xvec);
        static void CarttoSpherical (std::vector<std::vector<double>>& coords);
        static void SphericaltoCart (std::vector<std::vector<double>>& coords);
        
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
        
        static void InEFtoKS (std::vector<double>& xvec, double mass = 1,
                             const char* EFtype = "original");
        static void InEFtoKS (std::vector<std::vector<double>>& coords, 
                             double mass = 1,
                             const char* EFtype = "original");
        static void KStoInEF (std::vector<double>& xvec, double mass = 1,
                             const char* EFtype = "original");
        static void KStoInEF (std::vector<std::vector<double>>& coords,
                             double mass = 1,
                             const char* EFtype = "original");
        
        static void switchInEF (std::vector<double>& xvec, 
                                const char* from = "original");
        static void switchInEF (std::vector<std::vector<double>>& coords,
                                const char* from = "original");
        
        template <typename whatever1, typename whatever2>
        static void EF_from_uv_to_original
                    (std::vector<std::vector<double>>& coords, 
                     whatever1 a = 0, whatever2 b = 0)
                    {return Spacetime::switchInEF(coords, "uv");}
        
        template <typename whatever1, typename whatever2>
        static void EF_from_original_to_uv
                    (std::vector<std::vector<double>>& coords, 
                     whatever1 a = 0, whatever2 b = 0)
                    {return Spacetime::switchInEF(coords, "original");}
        
        template <typename whatever1, typename whatever2, typename whatever3>
        static void do_nothing(whatever1 a=0, whatever2 b=0, whatever3 c=0)
        {return;}

        // The General Causalities

        static std::vector<bool> BH_general_causal2D (const std::vector<double>& xvec, 
                                                const std::vector<double>& yvec,
                                                std::vector<double> period={},
                                                double mass = 1);

        static std::vector<bool> BH_general_causal3D (const std::vector<double>& xvec, 
                                                const std::vector<double>& yvec,
                                                std::vector<double> period={},
                                                double mass = 1); 

        static std::vector<bool> BH_general_causal4D (const std::vector<double>& xvec, 
                                                const std::vector<double>& yvec,
                                                std::vector<double> period={},
                                                double mass = 1);
        
        static std::vector<bool> BH_general_last_resort(const std::vector<double>& xvec, 
                                                const std::vector<double>& yvec,
                                                double mass = 1);
};

#endif /* SPACETIMES_H */