/// \authors Vid Homsak, Stefano Veroni
/// \date 28/09/2022

#ifndef SPRINKLEDCAUSET_H
#define SPRINKLEDCAUSET_H

#include <set>
#include <vector>

#include "causet.h"
#include "embeddedcauset.h"
#include "spacetimes.h"
#include "shapes.h"

using std::vector;
using std::set;

/**
 * @brief An causet obtained from sprinkling in a spacetime subset 
 *        of a specified shape.
 */
class SprinkledCauset: public EmbeddedCauset
{
    public:
        // Public Attributes
        double _intensity;
       
        // Constructor
        SprinkledCauset(int label = -1,
                    set<CausetEvent> past = {},
                    set<CausetEvent> future = {}, 
                    vector<double> position = {},
                    vector<double> coordinates = {});
        
        // Properties
        double Intensity();
        double Density(); //overwrites superclass?
        double LengthScale(); //overwrites superclass?

        // Methods
        vector<vector<double>> _sprinkle_coords();
        SprinkledCauset sprinkle();
        SprinkledCauset intensify();
        SprinkledCauset create();
        void add();
        void discard();      
};

#endif /* SPRINKLEDCAUSET_H */