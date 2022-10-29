/// \authors Vid Homsak, Stefano Veroni
/// \date 28/09/2022

#ifndef SPRINKLEDCAUSET_H
#define SPRINKLEDCAUSET_H

#include <set>
#include <vector>

#include "causetevent.h"
#include "causet.h"
#include "embeddedcauset.h"
#include "spacetimes.h"
#include "shapes.h"

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
                    std::set<CausetEvent> past = {},
                    std::set<CausetEvent> future = {}, 
                    std::vector<double> position = {},
                    std::vector<double> coordinates = {});
        
        // Properties
        double Intensity();
        double Density(); //overwrites superclass?
        double LengthScale(); //overwrites superclass?

        // Methods
        std::vector<std::vector<double>> _sprinkle_coords();
        std::set<CausetEvent> sprinkle();
        std::set<CausetEvent> intensify();
        std::set<CausetEvent> create();
        void add();
        void discard();
        
};

#endif /* SPRINKLEDCAUSET_H */