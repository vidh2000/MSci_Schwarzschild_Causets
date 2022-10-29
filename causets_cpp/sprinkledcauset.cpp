#ifndef SPRINKLEDCAUSET_H
#define SPRINKLEDCAUSET_H

#include <set>
#include <vector>

// For Vid (can't locate headers for some reason)
// Path: D:\Documents\Sola\Imperial College London\Year 4\MSci project\
            Project\causets_code\causets_cpp\"header".h...


#include <D:\Documents\Sola\Imperial College London\Year 4\MSci project\Project\causets_code\causets_cpp\spacetimes.h>


//////////////////////////////////////////////////////////////////////////////
// ::SprinkledCauset::
// Handles a causal set that is embedded in a subset of a manifold.
//////////////////////////////////////////////////////////////////////////////

SprinkledCauset::SprinkledCauset(int card,
                    double intensity,
                    int dim,
                    Spacetime spacetime,
                    std::vector<std::pair<std::string, CoordinateShape>> shape)
/**
 * @Constructor:
 * @brief Generates a sprinkled causal set by sprinkling in a spacetime subset. 
 *       The arguments `dim`, `shape` and `spacetime` are handled by the super 
 *       class `EmbeddedCauset` before events are sprinkled.
 *  Parameters
 *  -------------------
 *  @param card: int
 *      Number of sprinkled events.
 * @param intensity: float
 *      Sprinkling intensity parameter, the expected number of sprinkled 
 *      events.
 **/
{
    // initialise base class (EmbeddedCauset):
    //super().__init__(spacetime=spacetime, shape=shape, dim=dim)??????
    // sprinkle
    this -> intensity;
    if (card > 0)
    {sprinkle(card);}
    else
    {intensify(intensity);}
}

// Properties
double SprinkledCauset::Intensity()

    /* Returns the sprinkling intensity, which is the expected number of 
    sprinkled events. The exact number of sprinkled events is given by the 
    property 'Card'. */
    {return _intensity;}

double SprinkledCauset::Density()
    //overwrites superclass
    {return _intensity / *this.Shape().Volume();} 

double SprinkledCauset::LengthScale()
    //overwrites superclass
    {return std::pow(*this.Shape().Volume()/ _intensity, (1.0 / *this.Dim()));}

// Methods
std::vector<std::vector<double>> SprinkledCauset::_sprinkle_coords(
    int count, CoordinateShape shape)
{
    if (count < 0)
    {   throw std::invalid_argument(
            'The sprinkle cardinality has to be a non-negative integer.')
    }
    std::vector<std::vector<double>> coords;
    if ((shape.Name() == "cube") || (shape.Name() == "cuboid"))
    {
        std::vector<double> low;
        std::vector<double> high;
    
        low = shape.Center()-shape.Edges()/2;
        high = shape.Center()-shape.Edges()/2;

        // Generate coords randomly

        // Will be used to obtain a seed for the random number engine
        std::random_device rd;  
        // Standard mersenne_twister_engine seeded with rd()
        std::mt19937 gen(rd()); 
        std::uniform_real_distribution<> dis(low,high);
        for (int i=0; i<count;i++)
        {
            coords[i,:] = dis(gen);
        }
    }
    else if ((shape.Name() == "ball") || (shape.Name() == "cylinder") ||
             (shape.Name() == "diamond") || (shape.Name() == "bicone"))
    {
        // Create circle based sprinkle:
        bool isCylindrical = shape.Name()=="cylinder";
        bool isDiamond = (shape.Name()=="diamond") || (shape.Name()=="bicone")

        int d = *this.Dim()
        double b_r = shape.Parameter().radius;
        if (d==2 && isDiamond)
        {
            //pick "count" random coordinate tuples uniformly:
            std::vector<std::vector<double>> uv;
            // Random generator stuff
            std::random_device rd;  
            std::mt19937 gen(rd()); 
            std::uniform_real_distribution<> dis(-1.0,1.0);
            for (i=0;i<count;i++)
            {
                std::vector<double> = 
                for (j=0;j<2;j++)
                {
                    
                    uv[i]
                }
                
            } 
        }
    }
    
}


std::set<CausetEvent> SprinkledCauset::sprinkle(
                    int count, CoordinateShape shape=*this.Shape())
    /*
    Creates a fixed number of new events by sprinkling into `shape` (by 
        default the entire embedding region).
    */
{
    if (count<0)
    {   throw std::invalid_argument(
            'The sprinkle cardinality has to be a non-negative integer.')
    }
    _intensity += 1.0 * count;
    std::vector<std::vector<double>> coords = _sprinkle_coords(count,shape); //+rng parameter?
    return super().create(coords);
}   


std::set<CausetEvent> SprinkledCauset::intensify(
            double intensity, CoordinateShape shape = *this.Shape())
        /*
        Creates an expected number of new events by sprinkling into `shape` (by 
        default the entire embedding region). The expected number is determined 
        by the Poisson distribution with the given `intensity` parameter.
        */
{
    if (intensity < 0.0)
    {   throw std::invalid_argument(
            "The intensity parameter has to be a non-negative float.")
    }
    _intensity += intensity
    std::default_random_engine generator;
    std::poisson_distribution<int> distribution(intensity);
    int count = int count = distribution(generator);
    std::vector<std::vector<double>> coords = _sprinkle_coords(count,shape); // in python it's ...count shape,rng)!
    return super().create(coords); //?? incomplete
}


std::set<CausetEvent> SprinkledCauset::create( //?? incomplete
        Union[Iterable[List[float]],
        Iterable[np.ndarray],
        np.ndarray] coords,
        std::string labelformat = None, //(optional)
        bool relate = true)
{
    double card_old = 1.0* *this.Card();
    std::set<CausetEvent> eventSet = super().create(
                coords,labelFormat,relate);
    _intensity += *this.Card() * 1.0 - card_old;
    return eventSet;
}


void SprinkledCauset::add(std::set<CausetEvent> eventSet, bool unlink=false)
{
    double card_old = 1.0 * *this.Card();
    super().add(eventSet, unlink);  //??
    _intensity += (*this.Card() * 1.0 - card_old);
}


void SprinkledCauset::discard(std::set<CausetEvent> eventSet, bool unlink=false)
{
    double card_old = 1.0 * *this.Card();
    super().discard(eventSet, unlink);  //??
    _intensity *= (*this.Card() * 1.0 / card_old);
}
