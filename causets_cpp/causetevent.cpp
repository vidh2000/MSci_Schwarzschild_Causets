#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <fstream>
#include <stack>
#include <string>
#include <stdio.h>
#include <vector>
#include <set>
//#include <functions>

class CausetEvent
{
    /*
    Handles a single event (point) and its causal relations in a causal set.
    The attribute 'Label' can be used to assign a label, but does not need to 
    be unique (default: None).

    Instances of CausetEvent are comparable:
    a == b             is True if a is the same instance as b
    a < b              is True if a precedes b
    a <= b             is True if a precedes b or a is the same as b
    a > b              is True if a succeeds b
    a >= b             is True if a succeeds b or a is the same as b
    a.isSpacelikeTo(b) is True if a is spacelike separated to b
    
    */

    public:
    // Public Attributes
    std::string Label;
    std::set<CausetEvent> _prec;
    std::set<CausetEvent> _succ;
    std::vector<double> _coordinates;
    std::vector<double> _position;

    /** 
    * @Constructor:
    * Initialise a CausetEvent.
    * 
    *    Keyword parameters:
    * label: str
    *    Label for the event (does not need to be unique in a causet)
    * past: Iterable[CausetEvent]
    *    Set of past events (that may or may not be linked).
    *    This instance will automatically be added to their future.
    * future: Iterable[CausetEvent]
    *    Set of future events (that may or may not be linked).
    *    This instance will automatically be added to their past.
    * coordinates: Iterable[float]
    *    Coordinates if the event shall be considered as embedded in a 
    *    spacetime region.
    * position: Iterable[float]
    *    Coordinate pair of the event in a Hasse diagram if the Hasse 
    *    diagram is manually defined.
    **/
    CausetEvent(std::string label = "",
            std::set<CausetEvent> past = {},
            std::set<CausetEvent> future = {}, 
            std::vector<double> position = {},
            std::vector<double> coordinates = {}
            )
    {
    this -> Label = label;
    _prec = {}
    for (CausetEvent e = past.begin(); e != past.end(); e++){
        std::set<CausetEvent> presOrPast = e.PresentOrPast();
        _prec = set_add(_prec, presOrPast);} 
    _succ = {}  
    for (CausetEvent e = future.begin(); e != future.end();e++){
        std::set<CausetEvent> presOrFut = e.PresentOrFuture();
        _succ = set_add(_succ, presOrFut);}
    this -> _coordinates = coordinates;
    this -> _position = position;
    // Add this instance to its causal relatives:
    for (CausetEvent e = _prec.begin(); _prec.end(); e++)
    {   e._addToFuture(*this);}
    for (CausetEvent e = _succ.begin(); _succ.end(); e++)
    {   e._addToPast(*this);}

    
    //ACTIVE FUNCTIONALITIES
    //private:
    bool _addToPast(CausetEvent other){
        /*
        Adds an event to the past of this event.
        It returns False if the event is already in the past,
        otherwise it adds the event and return True.
        */
        const bool is_in = set_contains(other, _prec);
        if is_in
        {
            return false;
        }
        else
        {
            if (hasBeenLinked() && isLink(other, *this))
            {
                _lprec = set_diff(_lprec, other._lprec);
                _lprec.insert(other);
            }
            _prec.insert(other);
            return true;
        }
    }
    bool _addToFuture();
    bool _discard();

    // Other methods are public
    //public:

    //ACTIVE FUNCTIONALITIES
    void disjoin();
    void link();
    void unlink();
    // CausetEvent copy();
    // double Rank();

    // RELATIONAL BOOLEANS 
    bool hasBeenLinked();
    static bool isLink(); //@classmethod
    bool isPastLink();
    bool isFutureLink();
    bool isCausalTo();
    bool isLinkedTo();
    bool isSpacelikeTo();

    // PAST, FUTURE and SPACELIKE
    static int LinkCountOf(); //@staticmethod
    static std::set<CausetEvent> Past(); //@property
    static std::set<CausetEvent> Future(); //@property
    static std::set<CausetEvent> Cone(); //@property
    static std::set<CausetEvent> PresentOrPast(); //@property
    static std::set<CausetEvent> PresentOrFuture(); //@property
    std::set<CausetEvent> Spacelike();
    static int PastCard(); //@property
    static int FutureCard(); //@property
    static int ConeCard(); //@property
    int SpacelikeCard();
    static std::set<CausetEvent> LinkPast(); //@property
    static std::set<CausetEvent> LinkFuture(); //@property
    static std::set<CausetEvent> LinkCone(); //@property
    static int LinkPastCard(); //@property
    static int LinkFutureCard(); //@property
    static int LinkConeCard(); //@property

    // EMBEDDING FUNCTIONALITIES
    void embed();
    void disembed();
    static std::vector<double> Position(); //@property
    void SetPosition(); //@Poisition.setter
    static bool isEmbedded(); //@property
    static std::vector<double> Coordinates(); //@property
    static int CoordinatesDim(); //@property 

    
};









