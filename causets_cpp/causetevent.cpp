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

        // Public methods
        void disjoin;
        void link;
        void unlink;
        // CausetEvent copy;
        // double Rank;
        bool hasBeenLinked;
        
            //Class methods =  static??
            bool isLink;

        bool isPastLink;
        bool isFutureLink;
        bool isCausalTo;
        bool isLinkedTo;
        bool isSpacelikeTo;
        



    // Private Attributes
    std::set<CausetEvent> _prec;
    std::set<CausetEvent> _succ;
    std::vector<double> _coordinates;
    std::vector<double> _position;
    // Private methods
    bool _addToPast;
    bool _addToFuture;
    bool _discard;
    
    
    /*##### Constructor: ######
    Initialise a CausetEvent.

        Keyword parameters:
        label: str
            Label for the event (does not need to be unique in a causet)
        past: Iterable[CausetEvent]
            Set of past events (that may or may not be linked).
            This instance will automatically be added to their future.
        future: Iterable[CausetEvent]
            Set of future events (that may or may not be linked).
            This instance will automatically be added to their past.
        coordinates: Iterable[float]
            Coordinates if the event shall be considered as embedded in a 
            spacetime region.
        position: Iterable[float]
            Coordinate pair of the event in a Hasse diagram if the Hasse 
            diagram is manually defined.
    */
    CausetEvent(std::string label,
                std::set<CausetEvent> past,
                std::set<CausetEvent> future,
                std::vector<double> position;
                std::vector<double> coordinates)
    {
    
    this -> Label = label;
    //past
    //future
    this -> _coordinates = coordinates;
    this -> _position = position;
    // Add this instance to its causal relatives:
    /*for e in self._prec:
            e._addToFuture(self)
        for e in self._succ:
            e._addToPast(self)
    */
    } 

};











