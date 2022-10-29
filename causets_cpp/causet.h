/// \authors Vid Homsak, Stefano Veroni
/// \date 25/09/2022

#ifndef CAUSET_H
#define CAUSET_H

#include <cmath>
#include <set>
#include <vector>

#include "causetevent.h"

using namespace std;

class Causet
{
    public:
        set<CausetEvent> _events;

        // CONSTRUCTORS
        Causet(set<CausetEvent> set = {});
        Causet(vector<int> permutation); //FromPermutation     
        Causet(int n, 
                bool chain = false,           //NewChain
                bool antichain = false,       //NewAntiChain
                bool simplex = false,         //NewSimplex
                bool newcrown = false);  
        Causet(int length, int height, bool closed = true); //NewFence
        //NewKROrder not implemented
        Causet(vector<vector<double>>, bool fromCmatrix = true, 
                bool reverse_causality = false);
        
        // CAUSET MODIFICATIONS
        Causet merge(set<CausetEvent> pastSet, 
                     set<CausetEvent> futSet,
                     bool disjoint = false);
        Causet add(set<CausetEvent> eventSet, bool unlink = false);
        Causet discard(set<CausetEvent> eventSet, bool unlink = false);
        Causet coarsegrain(int card = nan(""), double perc = 0.2);
        Causet cgrain(int card = nan(""), double perc = 0.2);

        // CAUSET REPRESENTATION & SAVING
        vector<CausetEvent> nlabel(const char* method = "label", 
                                        bool reverse = false);
        vector<CausetEvent> nlist(const char* method = "label", 
                                        bool reverse = false);
        vector<vector<double>> CMatrix (const char* method="causality",
                                         set<CausetEvent> Events = {});
        vector<vector<double>> CMatrix (const char* method="causality",
                                         vector<CausetEvent> Events = {});
        vector<vector<double>> CTMatrix (const char* method="causality",
                                         set<CausetEvent> Events = {});
        vector<vector<double>> CTMatrix (const char* method="causality",
                                         vector<CausetEvent> Events = {});
        vector<vector<double>> LMatrix (const char* method="causality",
                                         set<CausetEvent> Events = {});
        vector<vector<double>> LMatrix (const char* method="causality",
                                         vector<CausetEvent> Events = {});
        vector<vector<double>> LTMatrix (const char* method="causality",
                                         set<CausetEvent> Events = {});
        vector<vector<double>> LTMatrix (const char* method="causality",
                                         vector<CausetEvent> Events = {});
        void saveCasCSV(const char* filename);
        void saveCasTXT(const char* filename);
        void saveLasCSV(const char* filename);
        void saveLasTXT(const char* filename);
        vector<CausetEvent> sortedByCausality(vector<CausetEvent> Events = {},
                                              bool reverse = false);
        vector<CausetEvent> sortedByLabel    (vector<CausetEvent> Events = {},
                                              bool reverse = false);

        // MATHEMATICS
        bool contains(CausetEvent x);
        Causet causet_union(Causet C);
        Causet causet_intersect(Causet C);
        Causet causet_diff(Causet C);
        Causet causet_symm_diff(Causet C);
        
        // KINEMATICS
        int size();
        int Card();
        double ord_fr(set<CausetEvent> A,
                      const char* denominator = "choose",
                      bool isdisjoined = true);
        double ord_fr(CausetEvent a, CausetEvent b,
                      const char* denominator = "choose",
                      bool isdisjoined = true);
        static double optimiser_placeholder();
        typedef double (*func)();//need to pick optimiser;
        double MMdim_est(const char* method = "random",
                        int d0 = 2, int Nsamples = 20,
                        int size_min = 10, double size_max = nan(""),
                        func optimiser = Causet::optimiser_placeholder);
        
        // LINKS
        // void link();
        // void unlink();
        // void LinkCount();

        // FIND
        CausetEvent find (int label);
        CausetEvent findAny (vector<int> labels);
        CausetEvent findAll (vector<int> labels);

        // PAST, PRESENT, FUTURE & NOTHING AT ALL
        static set<CausetEvent> PastOf (set<CausetEvent> Set, 
                                        bool includePresent = false,
                                        bool intersect = false);
        static set<CausetEvent> FutureOf (set<CausetEvent> Set, 
                                        bool includePresent = false,
                                        bool intersect = false);
        static set<CausetEvent> ConeOf (set<CausetEvent> Set, 
                                        bool includePresent = false,
                                        bool intersect = false);
        static set<CausetEvent> SpacelikeTo (set<CausetEvent> Set, 
                                        bool includePresent = false,
                                        bool intersect = false);
        
        // INTERVAL
        set<CausetEvent> Interval(CausetEvent a, CausetEvent b,
                                  bool includeBoundary = true,
                                  bool disjoin = false,
                                  const char* createmethod = "set");
        int IntervalCard(CausetEvent a, CausetEvent b,
                         bool includeBoundary = true);
        set<CausetEvent> PerimetralEvents(CausetEvent a, CausetEvent b,
                                        bool includeBoundary = true);
        int PerimetralCard(CausetEvent a, CausetEvent b,
                             bool includeBoundary = true);
        set<CausetEvent> InternalEvents(CausetEvent a, CausetEvent b,
                                        bool includeBoundary = true);
        int InternalCard(CausetEvent a, CausetEvent b,
                             bool includeBoundary = true);
        
        // CHAINS, ANTICHAINS & PATHS
        bool isChain(set<CausetEvent> Events);
        bool isAntiChain(set<CausetEvent> Events);
        set<set<CausetEvent>> Paths(CausetEvent a, CausetEvent b, 
                                    const char* length = "any");
        bool isPath(set<CausetEvent> Events);

        // PAST/FUTURE INFINITE
        set<CausetEvent> PastInf ();
        set<CausetEvent> FutureInf ();
        int PastInfCard ();
        int FutureInfCard ();
        set<CausetEvent> PastInfOf (set<CausetEvent> Events);
        set<CausetEvent> FutureInfOf (set<CausetEvent> Events);
        int PastInfCardOf (set<CausetEvent> Events);
        int FutureInfCardOf (set<CausetEvent> Events);
        set<CausetEvent> CentralAntichain (set<CausetEvent> Events);

        // OTHERS

};

#endif /*CAUSET_H*/