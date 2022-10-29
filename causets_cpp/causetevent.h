#ifndef CAUSETEVENT_H
#define CAUSETEVENT_H

#include <set>
#include <vector>

class CausetEvent
// Note no link fuctions have been included.
{
    public:

        // Public Attributes
        int Label;
        std::set<CausetEvent> _prec;
        std::set<CausetEvent> _succ;
        std::vector<double> _coordinates;
        std::vector<double> _position;

        // Constructor
        CausetEvent(int label = -1,
                    std::set<CausetEvent> past = {},
                    std::set<CausetEvent> future = {}, 
                    std::vector<double> position = {},
                    std::vector<double> coordinates = {});
        
        // "Getters"
        std::set<CausetEvent>& Past = _prec;
        std::set<CausetEvent>& Future = _succ;
        std::set<CausetEvent> Cone();
        std::set<CausetEvent> PresentOrPast();
        std::set<CausetEvent> PresentOrFuture();
        std::set<CausetEvent> Spacelike(std::set<CausetEvent> eventset);
        int PastCard();
        int FutureCard();
        int ConeCard();
        int SpacelikeCard(std::set<CausetEvent> eventset);
    
        //Overloading Operators
        bool operator ==(const CausetEvent &other) const;
        bool operator <(const CausetEvent &other) const;
        bool operator <=(const CausetEvent &other) const;
        bool operator >(const CausetEvent &other) const;
        bool operator >=(const CausetEvent &other) const;

        //Active Functionalities
        bool _addToPast(CausetEvent other);
        bool _addToFuture(CausetEvent other);
        bool _discard(CausetEvent other);
        void disjoin();

        //Relational Booleans
        bool isCausalTo(CausetEvent other);
        bool isSpacelikeTo(CausetEvent other);

        //Embedding Functionalities
        void embed(std::vector<double> coordinates, bool reembed = false);
        void disembed();
        bool isEmbedded();
        std::vector<double> Coordinates();
        int CoordinatesDim();

        // Hasse diagram functionalities
        std::vector<double> Position();
        void SetPosition(std::vector<double> value);





};

#endif /* CAUSETEVENT_H */