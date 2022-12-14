/// \authors Vid Homsak, Stefano Veroni
/// \date 29/10/2022

#ifndef EMBEDDEDCAUSET_H
#define EMBEDDEDCAUSET_H

#include <cmath>
#include <vector>
#include <unordered_set>

#include "causet.h"
#include "spacetimes.h"
#include "shapes.h"


/**
 * @brief An embedded causet in a spacetime subset of a specified shape.
 */
class EmbeddedCauset: public Causet
{
    public:
        //Inherited
        // vector<vector<int8_t>> _CMatrix = {};
        // bool _special_matrix = false;
        // int _size = 0;
        // int _dim = 0;
        
        // vector<std::unordered_set<int>> _pasts   = {};
        // vector<std::unordered_set<int>> _futures = {};
        // vector<std::unordered_set<int>> _past_links   = {};
        // vector<std::unordered_set<int>> _future_links = {};
        std::vector<std::vector<double>> _coords = {};
        CoordinateShape _shape = CoordinateShape();
        Spacetime _spacetime = Spacetime();

        // CONSTRUCTOR
        
        EmbeddedCauset();
        EmbeddedCauset(Spacetime spacetime, 
                       CoordinateShape shape, 
                       std::vector<std::vector<double>> coordinates,
                       bool make_matrix = true,
                       bool special = false,
                       bool use_transitivity = true,
                       bool make_sets = false,
                       bool make_links = true,
                       const char* generation_mode = "future");
        

        // GETTERS

        int spacetime_dim();
        double density();
        double length_scale();
        CoordinateShape &get_shape = _shape;
        Spacetime &get_spacetime = _spacetime;

        std::vector<double> eu_distances();
        std::vector<double> sp_radii();
        double max_eu_dist();
        double max_sp_rad();
        double max_along(int ndim = 0);
        double min_along(int ndim = 0);


        // CAUSAL RELATIONS

        bool causality(std::vector<double> xvec, 
                       std::vector<double> yvec);
        std::vector<bool> general_causality(std::vector<double> xvec, 
                                            std::vector<double> yvec);
        bool areTimelike4D(std::vector<double> &xvec, std::vector<double> &yvec);
        bool AprecB(std::vector<double> xvec, std::vector<double> yvec);


        // MODIFIERS (BETA: NONE TESTED)

        void sort_coords(int dim = 0, bool reverse = false);
        void relabel(const char* method = "0", bool reverse = false);   
        void add(std::vector<double> xvec);
        void discard(int label, bool make_matrix = true, 
                     bool make_sets = false, bool make_links = true);  
        void discard(std::vector<int> labels, bool make_matrix = true, 
                     bool make_sets = false, bool make_links = true);

        // SAVE
        
        void save_causet(const char* path_file_ext,
                         const char* storage_option = "sets");


        //MAKE ATTRIBUTES
        
        void make_attrs(const char* method = "coordinates",
                            bool make_matrix = true,
                            bool special = false,
                            bool use_transitivity = true,
                            bool make_sets = false,
                            bool make_links = true,
                            const char* sets_type = "future");

        // Behind the scenes

        void make_all_but_links();
        void make_cmatrix(const char* method = "coordinates",
                            bool special = true,
                            bool use_transitivity = true);
        void make_cmatrix_and_allpasts(bool special = true);
        void make_cmatrix_and_allfuts(bool special = true);
        void make_cmatrix_and_pasts(const char* method = "coordinates",
                                        bool special = true,
                                        bool use_transitivity = true);
        void make_cmatrix_and_futs(const char* method = "coordinates",
                                        bool special = true,
                                        bool use_transitivity = true);
        void make_cmatrix_and_pastlinks(const char* method = "coordinates",
                                        bool special = true);
        void make_cmatrix_and_futlinks(const char* method = "coordinates",
                                        bool special = true);
        void make_sets          (const char* method = "coordinates");
        void make_all_pasts     (const char* method = "coordinates");
        void make_all_futures   (const char* method = "coordinates");
        void make_pasts         (const char* method = "coordinates");
        void make_futures       (const char* method = "coordinates");
        void make_past_links    (const char* method = "coordinates");
        void make_fut_links     (const char* method = "coordinates");
        
        // Methods for counting links -> entropy
        
        int count_links_fromCMatrix(double& t_f, double r_S = 2);
        int count_links_BH(double& t_f, double r_S = 2);

        std::map<int,std::vector<int>> get_lambdas_fromCMatrix(double& t_f, 
                                                               double r_S = 2);
        std::map<int,int> count_lambdas_fromCMatrix(double& t_f, double r_S = 2);
        void save_lambdas(const char* path_file_ext, const char* storage_option, 
                            double & t_f, double & r_S);
        

        // Counting Behind the Scenes

        std::map<int,std::vector<int>> get_lambdas_from_futlinks(double& t_f,double r_S = 2);
        std::map<int,int> get_lambdas_sizes(double& t_f, double r_S = 2);
        std::map<int,int> get_lambdas_distr(const std::map<int, int> & lambdas);


        //Destructor
        ~EmbeddedCauset();       
};

#endif /*EMBEDDEDCAUSET_H*/