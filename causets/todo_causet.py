#!/usr/bin/env python
'''
Created on 20 Jul 2020

@author: Christoph Minz
@license: BSD 3-Clause

Stefano Veroni added some new functions 
(those with """ """ description rather than ')
'''
from __future__ import annotations
from typing import Set, Iterable, List, Any, Tuple, Iterator, Union, Optional
from causets.causetevent import CausetEvent  # @UnresolvedImport

import numpy as np
import random
import itertools


class Causet(object):
    '''
    Causal set class to handle operations of a set of `CausetEvent`.
    '''

    _events: Set[CausetEvent]

    def __init__(self, eventSet: Set[CausetEvent] = set()) -> None:
        '''
        Generates a Causet class instance from a set of `CausetEvent`. The 
        `CausetEvent` instances are not checked for logical consistency.
        '''
        while True:
            l: int = len(eventSet)
            eventSet = Causet.ConeOf(eventSet)
            if len(eventSet) == l:
                break
        self._events: Set[CausetEvent] = eventSet

    def __iter__(self) -> Iterator[CausetEvent]:
        return iter(self._events)

    def __repr__(self) -> str:
        return repr(self._events)


    ####################################################################
    #___________________________________________________________________
    # "CREATORS" OF CAUSETS
    #___________________________________________________________________
    ####################################################################


    @staticmethod
    def FromCausalMatrix(C: np.ndarray) -> 'Causet':
        """
        Converts a logical matrix into a `Causet` object. The entry `C[i, j]` 
        has to be True or 1 if the event with index j is in the (link) past of 
        event with index i. If the matrix has less rows than columns, empty 
        rows are added after the last row. However, if the matrix has more rows 
        than columns, a ValueError is raised. A ValueError is also raised if 
        the matrix contains causal loops.
        """
        rowcount: int = C.shape[0]
        colcount: int = C.shape[1]
        if colcount < rowcount:
            raise ValueError('The specified matrix cannot be extended ' +
                             'to a square matrix.')
        events: np.ndarray = np.array([CausetEvent(label=i)
                                       for i in range(1, colcount + 1)])
        e: CausetEvent
        for i in range(rowcount):
            e = events[i]
            future: Set[CausetEvent] = set(events[np.where(C[i, :])[0]])
            past: Set[CausetEvent] = set(events[np.where(C[:, i])[0]])
            if (past & future) or (e in past) or (e in future):
                raise ValueError('The causet is not anti-symmetric.')
            e._prec = past
            e._succ = future
        # complete pasts and futures (if the input contains links only)
        for i in range(colcount):
            e = events[i]
            e._prec = Causet.PastOf(e._prec, includePresent=True)
            e._succ = Causet.FutureOf(e._succ, includePresent=True)
        return Causet(set(events))

    @staticmethod
    def FromTCausalMatrix(C: np.ndarray) -> 'Causet':
        """
        Returns `FromPastMatrix` of the transposed input. 
        """
        return Causet.FromPastMatrix(C)

    @staticmethod #replaced by FromCausalMatrix
    def FromPastMatrix(C: np.ndarray) -> 'Causet':
        '''
        Converts a logical matrix into a `Causet` object. The entry `C[i, j]` 
        has to be True or 1 if the event with index j is in the (link) past of 
        event with index i. If the matrix has less rows than columns, empty 
        rows are added after the last row. However, if the matrix has more rows 
        than columns, a ValueError is raised. A ValueError is also raised if 
        the matrix contains causal loops.
        '''
        rowcount: int = C.shape[0]
        colcount: int = C.shape[1]
        if colcount < rowcount:
            raise ValueError('The specified matrix cannot be extended ' +
                             'to a square matrix.')
        events: np.ndarray = np.array([CausetEvent(label=i)
                                       for i in range(1, colcount + 1)])
        e: CausetEvent
        for i in range(rowcount):
            e = events[i]
            past: Set[CausetEvent] = set(events[np.where(C[i, :])[0]])
            future: Set[CausetEvent] = set(events[np.where(C[:, i])[0]])
            if (past & future) or (e in past) or (e in future):
                raise ValueError('The causet is not anti-symmetric.')
            e._prec = past
            e._succ = future
        # complete pasts and futures (if the input contains links only)
        for i in range(colcount):
            e = events[i]
            e._prec = Causet.PastOf(e._prec, includePresent=True)
            e._succ = Causet.FutureOf(e._succ, includePresent=True)
        return Causet(set(events))

    @staticmethod #replaced AND IMPROVED by FromTCausalMatrix
    def FromFutureMatrix(C: np.ndarray) -> 'Causet':
        '''
        Returns `FromPastMatrix` of the transposed input. 
        '''
        return Causet.FromPastMatrix(C.T)



    ####################################################################
    #___________________________________________________________________
    # SET MODIFICATION
    #___________________________________________________________________
    ####################################################################
    @staticmethod
    def merge(pastSet: Iterable, futureSet: Iterable,
              disjoint: bool = False) -> 'Causet':
        '''
        Returns a new Causet instance that joins the event sets `pastSet` and 
        `futureSet`.
        If not `disjoint` (default), then the event of `pastSet` are also 
        assigned to the the past of every event in `futureSet` and vice versa.
        '''
        if not disjoint:  # add pastSet as past of futureSet
            for p in pastSet:
                for f in futureSet:
                    p._addToFuture(f)
                    f._addToPast(p)
        return Causet(set(pastSet) | set(futureSet))

    def add(self, eventSet: Iterable, unlink: bool = False) -> None:
        '''
        Adds all the events of the (causal) set `eventSet` (`Causet` or 
        `Set[CausetEvent]`) to this instance. 
        '''
        self._events.update(eventSet)
        if unlink:
            for e in self._events:
                e.unlink()
        if hasattr(self, '__diagram_coords'):
            delattr(self, '__diagram_coords')

    def discard(self, eventSet: Iterable, unlink: bool = False) -> None:
        '''
        Discards all the events of the (causal) set `eventSet` (`Causet` or 
        `Set[CausetEvent]`) from this instance. 
        '''
        self._events.difference_update(eventSet)
        if unlink:
            for e in self._events:
                e.unlink()
        if hasattr(self, '__diagram_coords'):
            delattr(self, '__diagram_coords')

    def coarsegrain(self, card = None, perc : float = 0.2):
        """
        Coarse Grain the Causet by:
        - removing <card> events, if card given;
        - removing <int(perc * C.Card)> events, if card not given;\n
        Default is card None and perc = 0.2.
        """
        if card is None:
            card = int(len(self) * perc)
        events = random.sample(self._events, k = card)
        self.discard(events)
    
    def cgrain(self, card = None, perc : float = 0.2):
        self.coarsegrain(card = card, perc = perc)



    ####################################################################
    #___________________________________________________________________
    # CAUSET REPRESENTATION & SAVING TO FILE
    #___________________________________________________________________
    ####################################################################
    def nlabel (self, method = "label", reverse = False) -> List[CausetEvent]:
        """
        Turn Causet into list of events with a natural labelling. Methods are:
        - 'label': use labels (if created from causal matrix definitely best as
        it preserves the indexing of the matrix, though note label is 1 not 0)
        - 'past': by listing the set and sorting by cardinality of past
        - 'causality': use SortedByCausality
        """
        a = 1 if not reverse else -1      
        if method == "label":
            Clist = list(self)
            Clist.sort(key = lambda e : a*e.Label)
        elif method == "label":
            Clist = list(self)
            Clist.sort(key = lambda e : a*e.PastCard)
        elif method == "causality":
            Clist = list(self)
            Clist = self.sortedByCausality(reverse)
        return Clist
    
    def nlist (self, method = "label", reverse = False) -> List[CausetEvent]:
        """
        Turn Causet into list of events with a natural labelling. Methods are:
        - 'label': use labels (if created from causal matrix definitely best as
        it preserves the indexing of the matrix, though note label is 1 not 0)
        - 'past':      by listing the set and sorting by cardinality of past
        - 'causality': use SortedByCausality
        """
        a = 1 if not reverse else -1      
        if method == "label":
            Clist = list(self)
            Clist.sort(key = lambda e : a*e.Label)
        elif method == "label":
            Clist = list(self)
            Clist.sort(key = lambda e : a*e.PastCard)
        elif method == "causality":
            Clist = list(self)
            Clist = self.sortedByCausality(reverse)
        return Clist

    def CMatrix(self, 
                labeledEvents: Optional[List[CausetEvent]] = None,
                dtype: Any = int,
                method = "causality") -> np.ndarray:
        """
        Returns the logical causal matrix such that `C[i, j]` is 1 if 
        i preceds. 
        Events are indexed by `labeledEvents` (by default sorted by 
        causality).
        """
        if labeledEvents is None:
            labeledEvents = self.nlist(method)
        l: int = len(labeledEvents)
        C: np.ndarray = np.zeros((l, l), dtype)
        for i, a in enumerate(labeledEvents):
            for j, b in enumerate(labeledEvents):
                C[i, j] = 1 if a < b else 0
        return C
    
    def CTMatrix(self, 
                labeledEvents: Optional[List[CausetEvent]] = None,
                dtype: Any = int) -> np.ndarray:
        """
        Returns the logical causal PAST matrix such that `C[i, j]` is 1 if 
        j preceds i. 
        The events are indexed by `labeledEvents` (by default sorted by 
        causality).
        """
        if labeledEvents is None:
            labeledEvents = self.nlist()
        l: int = len(labeledEvents)
        C: np.ndarray = np.zeros((l, l), dtype)
        for i, a in enumerate(labeledEvents):
            for j, b in enumerate(labeledEvents):
                C[i, j] = 1 if a > b else 0
        return C

    
    # ENCAPSULATED INTO NLIST, BUT STILL VERY NECESSARY
    def sortedByLabels(self, eventSet: Optional[Set[CausetEvent]] = None,
                       reverse: bool = False) -> List[CausetEvent]:
        '''
        Returns the causet events `eventSet` (if None than the entire 
        causet) as a list sorted ascending (default) or descending by 
        their labels.
        '''
        if eventSet is None:
            eventSet = self._events
        unsortedList: List[CausetEvent] = list(eventSet)
        sortedIndex: np.ndarray = \
            np.argsort([e.Label for e in unsortedList])
        sortedList: List[CausetEvent] = \
            [unsortedList[i] for i in sortedIndex]
        if reverse:
            sortedList.reverse()
        return sortedList

    def sortedByCausality(self, eventSet: Optional[Set[CausetEvent]] = None,
                          reverse: bool = False) -> List[CausetEvent]:
        '''
        Returns the causet events `eventSet` (if None than the entire causet) 
        as a list sorted ascending (default) or descending by their causal 
        relations.
        '''
        if eventSet is None:
            eventSet = self._events
        eList: List[CausetEvent] = self.sortedByLabels(eventSet, reverse=True)
        c: int = len(eList)
        for i in range(c):
            for j in range(i):
                if eList[i] < eList[j]:
                    eList[i], eList[j] = eList[j], eList[i]
        if reverse:
            eList.reverse()
        return eList    


   
    ####################################################################
    #___________________________________________________________________
    # PROPER CAUSAL SETS KINEMATICS ANALYSIS
    #___________________________________________________________________
    ####################################################################
    def __len__(self):
        """
        len method returning cardinality
        """
        return len(self._events)

    @staticmethod
    def len(other: 'Causet') -> int:
        '''
        Returns the number of events (set cardinality) of some Causet instance. 
        '''
        return len(other._events)

    @property
    def Card(self) -> int:
        '''
        Returns the number of events (set cardinality) in this instance.
        '''
        return len(self._events)

    def ord_fr_A(self, A: Set[CausetEvent], den = "choose",
                isdisjoined = True) -> float:
        """
        Find the ordering fraction of an Aexandrov Interval: 
        the ratio of actual relations over possible in such interval.

        Parameters:
        -----------
        - A: list of [CausetEvent]\n
            Alexandrov Interval of which computing the ordering fraction

        - den: string ('choose' or 'n2')
            Use as denominator:
            - 'choose' -> |A|(|A|-1)/2, i.e. |A| choose 2 (Default).
            - 'n2'     -> (|A|^2)/2.
        
        - isdisjoined: Bool \n
            If True, A is taken as having been disjoined from rest of
            causet, meaning forall e in A e.Causal contains only events
            in A.

        Return:
        .......
        - Ordering fraction of Alexandrov Interval\n
            This is nrelations / (N choose 2)
        """
        if den != 'choose' and den != 'n2':
            raise ValueError ("den value is neither 'choose' nor 'n2'")
        N = len(A)
        nrelations = 0
        if isdisjoined:
            for ei in A:
                nrelations += ei.PastCard
        else:
            for ei in A:
                nrelations += len(ei.Past & A)
        fr = 2 * nrelations / ( N * (N - (den=='choose')) )
        return fr

    def ord_fr_ab(self, a: CausetEvent, b: CausetEvent,
                  den = "choose") -> float:
        """
        Find the ordering fraction of an interval between a and b: 
        the ratio of actual relations over possible in such interval.

        Parameters:
        -----------
        - a: CausetEvent\n
        
        - b: CausetEvent\n

        - mode: string
            Use as denominator:
            - 'choose' -> |A|(|A|-1)/2, i.e. |A| choose 2 (Default).
            - 'n2'     -> (|A|^2)/2.

        Return:
        .......
        - Ordering fraction of Alexandrov Interval
            This is nrelations / (N choose 2)
        """
        if den != 'choose' and den != 'n2':
            raise ValueError ("den value is neither 'choose' nor 'n2'")
        afutr = a.PresentOrFuture
        bpast = b.PresentOrPast
        A = afutr & bpast
        N = len(A)
        nrelations = 0
        for ei in A:
            nrelations += len(ei.Future & A)
        fr = 2 * nrelations / ( N * (N - (den=='choose')) )
        return fr


    ####################################################################
    #___________________________________________________________________
    # SET MATHEMATICS
    #___________________________________________________________________
    ####################################################################
    def __contains__(self, other: CausetEvent) -> bool:
        return other in self._events

    def __sub__(self, other: Iterable[CausetEvent]) -> Set[CausetEvent]:
        return self._events - set(other)

    def __rsub__(self, other: Iterable[CausetEvent]) -> Set[CausetEvent]:
        return self._events - set(other)

    def __or__(self, other: Iterable[CausetEvent]) -> Set[CausetEvent]:
        return self._events | set(other)

    def __ror__(self, other: Iterable[CausetEvent]) -> Set[CausetEvent]:
        return self._events | set(other)

    def __and__(self, other: Iterable[CausetEvent]) -> Set[CausetEvent]:
        return self._events & set(other)

    def __rand__(self, other: Iterable[CausetEvent]) -> Set[CausetEvent]:
        return self._events & set(other)

    def __xor__(self, other: Iterable[CausetEvent]) -> Set[CausetEvent]:
        return self._events ^ set(other)

    def __rxor__(self, other: Iterable[CausetEvent]) -> Set[CausetEvent]:
        return self._events ^ set(other)

    def difference(self, other: Iterable[CausetEvent]) -> Set[CausetEvent]:
        return self._events.difference(set(other))

    def intersection(self, other: Iterable[CausetEvent]) -> Set[CausetEvent]:
        return self._events.intersection(set(other))

    def symmetric_difference(self, other: Iterable[CausetEvent]) -> \
            Set[CausetEvent]:
        return self._events.symmetric_difference(set(other))

    def union(self, other: Iterable[CausetEvent]) -> Set[CausetEvent]:
        return self._events.union(set(other))
    


    ###################################################################
    #__________________________________________________________________
    # PAST, PRESENT, FUTURE & NOTHING AT ALL
    #__________________________________________________________________
    ###################################################################
    @staticmethod
    def PastOf(eventSet: Set[CausetEvent], includePresent: bool = False,
               intersect: bool = False) -> Set[CausetEvent]:
        '''
        Returns the set of events that are in the past of `eventSet`.
        '''
        newEventSet: Set[CausetEvent] = set()
        if includePresent and intersect:
            for e in eventSet:
                newEventSet &= e.PresentOrPast
        elif intersect:
            for e in eventSet:
                newEventSet &= e.Past
        else:
            for e in eventSet:
                newEventSet |= e.Past
            if includePresent:
                newEventSet |= eventSet
        return newEventSet

    @staticmethod
    def FutureOf(eventSet: Set[CausetEvent], includePresent: bool = False,
                 intersect: bool = False) -> Set[CausetEvent]:
        '''
        Returns the set of events that are in the future of `eventSet`.
        '''
        newEventSet: Set[CausetEvent] = set()
        if includePresent and intersect:
            for e in eventSet:
                newEventSet &= e.PresentOrFuture
        elif intersect:
            for e in eventSet:
                newEventSet &= e.Future
        else:
            for e in eventSet:
                newEventSet |= e.Future
            if includePresent:
                newEventSet |= eventSet
        return newEventSet

    @staticmethod
    def ConeOf(eventSet: Set[CausetEvent], includePresent: bool = True,
               intersect: bool = False) -> Set[CausetEvent]:
        '''
        Returns the set of events that are in the cone of `eventSet`.
        '''
        newEventSet: Set[CausetEvent] = set()
        if includePresent and intersect:
            for e in eventSet:
                newEventSet &= e.Cone
        elif intersect:
            for e in eventSet:
                newEventSet &= (e.Past | e.Future)
        else:
            for e in eventSet:
                newEventSet |= e.Past | e.Future
            if includePresent:
                newEventSet |= eventSet
        return newEventSet

    def SpacelikeTo(self, eventSet: Set[CausetEvent]) -> Set[CausetEvent]:
        '''
        Returns the set of events that are spacelike separated to `eventSet`.
        '''
        return self._events - self.ConeOf(eventSet, includePresent=True)


    ###################################################################
    #__________________________________________________________________
    # INTERVALS
    #__________________________________________________________________
    ################################################################### 
    @staticmethod
    def Interval(a: CausetEvent, b: CausetEvent,
                includeBoundary: bool = True,
                disjoin: bool = False,
                createmethod = "set") -> Set[CausetEvent]:
        """
        Returns the causal interval (Alexandrov set) between events 
        `a` and `b` or an empty set if not `a <= b`.

        Parameters
        -----------
        - a: CausetEvent\n
        
        - b: CausetEvent\n

        - includeBoundary: Bool\n
            If True (default), the events `a` and `b` are 
            included in the interval.

        - disjoin: Bool\n
            If True (non default), events lose information about
            the rest of the causet (and also embedding information).
        
        - createmethod: str\n
            - 'set': create from intersection of a.PresentOrFuture,\
                    b.PresentOrPast and the Causet. Default. Only option\
                    for disjoin == False. \n
            - 'matrix': create by creating and slicing CausalMatrix.\n
        """
        if createmethod != "set" and createmethod != "matrix":
            raise AttributeError("createmethod neither 'set' nor 'matrix")

        if not a <= b:
            AinC = set()
        elif a == b:
            AinC = {a}
        elif includeBoundary:
            AinC = a.PresentOrFuture & b.PresentOrPast
        else:
            AinC = a.Future & b.Past
        
        if disjoin and createmethod == "set":
            A = set()
            for einC in AinC:
                l = einC.Label 
                if l is None:
                    e = CausetEvent()
                elif isinstance(l, int):
                    e = CausetEvent(label = l)
                else:
                    e = CausetEvent(label = f'{l}')
                e._prec = einC._prec & a.PresentOrFuture
                e._succ = einC._succ & b.PresentOrPast
                A.add(e)

        elif disjoin and createmethod == "matrix":
            Am = AinC.CMatrix('label')\
                [int(a.Label):int(b.Label)+1, int(a.Label):int(b.Label)+1]
            useless_es = []
            for k in range(1, len(Am[0])-1):
                if Am[0][k]!=1 or Am[k][-1]!=1:
                    useless_es.append(k)
            if not includeBoundary:
                useless_es.append(0, len(Am[0]))
            Am = np.delete(Am, useless_es, 0)
            Am = np.delete(Am, useless_es, 1)
            A = Causet().FromCausalMatrix(Am)._events
        else:
            A = AinC
        return A
    
    @staticmethod
    def IntervalCard(a: CausetEvent, b: CausetEvent,
                     includeBoundary: bool = True) -> int:
        '''
        Returns the cardinality of the causal interval (Alexandrov set) between 
        events `a` and `b` or 0 if not `a <= b`.
        If `includeBoundary == True` (default), the events `a` and `b` are 
        included in the interval.
        '''
        if not a <= b:
            return 0
        elif a == b:
            return 1
        else:
            return len(a.Future & b.Past) + 2 * int(includeBoundary)
