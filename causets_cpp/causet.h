/// \authors Vid Homsak, Stefano Veroni
/// \date 25/09/2022

#ifndef CAUSET_H
#define CAUSET_H

#include <cmath>
#include <set>
#include <vector>

using std::vector;
using std::set;

class Causet
{
    public:
        vector<set<int>> _future;
        vector<set<int>> _past;

        // CONSTRUCTORS
        Causet(vector<set<int>> future = {}, 
                vector<set<int>> past = {});
        Causet(vector<vector<int>> M, bool Lmatrix = true, 
                bool reverse_causality = false);
        Causet(vector<vector<bool>> M, bool Lmatrix = true, 
                bool reverse_causality = false);
        // The following are included but not implemented
        Causet(vector<int> permutation); //FromPermutation     
        Causet(int n, 
                bool chain = false,           //NewChain
                bool antichain = false,       //NewAntiChain
                bool simplex = false,         //NewSimplex
                bool newcrown = false);  
        Causet(int length, int height, bool closed = true); //NewFence
        //NewKROrder not implemented
        
        // CAUSET MODIFICATIONS
        Causet merge(set<int> pastSet, 
                     set<int> futSet,
                     bool disjoint = false);
        Causet add(set<int> eventSet, bool unlink = false);
        Causet discard(set<int> eventSet, bool unlink = false);
        Causet coarsegrain(int card = nan(""), double perc = 0.2);
        Causet cgrain(int card = nan(""), double perc = 0.2);

        // CAUSET REPRESENTATION & SAVING
        vector<vector<double>> CMatrix (const char* method="causality",
                                         set<int> labels = {});
        vector<vector<double>> CMatrix (const char* method="causality",
                                         set<int> labels = {});
        vector<vector<double>> CTMatrix (const char* method="causality",
                                         set<int> labels = {});
        vector<vector<double>> CTMatrix (const char* method="causality",
                                         set<int> labels = {});
        vector<vector<double>> LMatrix (const char* method="causality",
                                         set<int> labels = {});
        vector<vector<double>> LMatrix (const char* method="causality",
                                         set<int> labels = {});
        vector<vector<double>> LTMatrix (const char* method="causality",
                                         set<int> labels = {});
        vector<vector<double>> LTMatrix (const char* method="causality",
                                         set<int> labels = {});
        void saveCasCSV(const char* filename);
        void saveCasTXT(const char* filename);
        void saveLasCSV(const char* filename);
        void saveLasTXT(const char* filename);
        vector<vector<set<int>>> sortedByCausality(set<int> labels = {},
                                                   bool reverse = false);
        vector<vector<set<int>>> sortedByCausality(set<int> labels = {},
                                                   bool reverse = false);

        
        // CAUSET EVENTS METHODS
        bool prec(int a, int b);
        bool spacelike(int a, int b);
        bool islink(int a, int b);


        // MATHEMATICS
        bool contains(int label = 0);
        bool contains(set<int> labels = {0});
        Causet causet_union(Causet C);
        Causet causet_intersect(Causet C);
        Causet causet_diff(Causet C);
        Causet causet_symm_diff(Causet C);

        
        // KINEMATICS
        int size();
        int Card();
        double ord_fr(Causet A,
                      const char* denominator = "choose",
                      bool isdisjoined = true);
        double ord_fr(vector<set<int>> A_future = {}, 
                      vector<set<int>> A_past = {},
                      const char* denominator = "choose",
                      bool isdisjoined = true);
        double ord_fr(int a, int b,
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
        // int find (int label);
        // int findAny (vector<int> labels);
        // int findAll (vector<int> labels);


        // OVERALL PAST, PRESENT, FUTURE & NOTHING AT ALL
        set<int> PastOf (int label, 
                                bool includePresent = false,
                                bool intersect = false);
        set<int> PastOf (set<int> labels, 
                                bool includePresent = false,
                                bool intersect = false);
        set<int> FutureOf (int label, 
                                  bool includePresent = false,
                                  bool intersect = false);
        set<int> FutureOf (set<int> labels, 
                                  bool includePresent = false,
                                  bool intersect = false);
        set<int> ConeOf (int label, 
                                bool includePresent = false,
                                bool intersect = false);
        set<int> ConeOf (set<int> labels, 
                                bool includePresent = false,
                                bool intersect = false);
        set<int> SpacelikeTo (int label, 
                                     bool includePresent = false,
                                     bool intersect = false);
        set<int> SpacelikeTo (set<int> labels, 
                                     bool includePresent = false,
                                     bool intersect = false);
        

        // INTERVAL
        Causet Interval(int a, int b,
                        bool includeBoundary = true,
                        bool disjoin = false,
                        const char* createmethod = "set");
        int IntervalCard(int a, int b, bool includeBoundary = true);
        Causet PerimetralEvents(int a, int b,  bool includeBoundary = true);
        int PerimetralCard(int a, int b, bool includeBoundary = true);
        Causet InternalEvents(int a, int b, bool includeBoundary = true);
        int InternalCard(int a, int b, bool includeBoundary = true);
        

        // CHAINS, ANTICHAINS & PATHS
        bool isChain(set<int> Events);
        bool isAntiChain(set<int> Events);
        set<set<int>> Paths(int a, int b, 
                            const char* length = "any");
        bool isPath(set<int> Events);


        // PAST/FUTURE INFINITE
        set<int> PastInf ();
        set<int> FutureInf ();
        int PastInfCard ();
        int FutureInfCard ();
        set<int> PastInfOf (set<int> Events = {});
        set<int> FutureInfOf (set<int> Events = {});
        int PastInfCardOf (set<int> Events = {});
        int FutureInfCardOf (set<int> Events = {});
        set<int> CentralAntichain (set<int> Events = {});

        // OTHERS

};

#endif /*CAUSET_H*/