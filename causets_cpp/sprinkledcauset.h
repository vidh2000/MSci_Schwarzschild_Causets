#ifndef SPRINKLEDCAUSET_H
#define SPRINKLEDCAUSET_H

#include <set>
#include <vector>

// For Vid (can't locate headers for some reason)
// Path: D:\Documents\Sola\Imperial College London\Year 4\MSci project\Project\causets_code\causets_cpp\"header".h...

#include <D:\Documents\Sola\Imperial College London\Year 4\MSci project\Project\causets_code\causets_cpp\spacetimes.h>

class SprinkledCauset
// Note no link fuctions have been included.
{
    public:
        // Public Attributes
        double _intensity;
       
        // Constructor
        SprinkledCauset(int card = 0,
                    double intensity = 0.0,
                    int dim = -1,
                    Spacetime spacetime,
                    std::vector<std::pair<std::string, CoordinateShape>> shape
                    );
        
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