#!/usr/bin/env python
'''
Created on 20 Jul 2020

@author: Christoph Minz
@license: BSD 3-Clause

Stefano Veroni added some new functions:
- nlabel and nlist
- ptime, MMdim_est
(those with """ """ description rather than ')
'''
from __future__ import annotations
from typing import Set, List, Iterable, Tuple  # @UnusedImport
from typing import Callable, Union, Optional  # @UnusedImport

import numpy as np
import scipy as sp
from scipy.optimize import fsolve
import random

from causets.causetevent import CausetEvent  # @UnresolvedImport
from causets.causet import Causet  # @UnresolvedImport
from causets.shapes import CoordinateShape  # @UnresolvedImport
from causets import spacetimes  # @UnresolvedImport
from causets.spacetimes import Spacetime  # @UnresolvedImport


class EmbeddedCauset(Causet):
    '''
    Handles a causal set that is embedded in a subset of a spacetime manifold.
    '''

    _shape: CoordinateShape
    _spacetime: Spacetime

    @staticmethod
    def __raiseDimValueError__(argument: str):
        e: str = 'EmbeddedCauset: The dimension of `%s` is ' + \
                 'not compatible with the other arguments.'
        raise ValueError(e % argument)

    def __init__(self,
                 spacetime: Optional[Union[Spacetime, str]] = None,
                 shape: Optional[Union[str, CoordinateShape]] = None,
                 coordinates: Optional[Union[List[List[float]],
                                             List[np.ndarray],
                                             np.ndarray]] = None,
                 dim: int = -1) -> None:
        '''
        Generates an embedded causal set in a spacetime subset of a specified 
        (coordinate) shape and with the events specified by `causet` and 
        `coordinates`.

        Optional parameters
        -------------------
        spacetime: Spacetime, str
            A spacetime object (including parameters) that determines the 
            causality, or name of spacetime to initialise with default 
            parameters. Supported values for the name are 'flat', 'Minkowski', 
            'dS', 'de Sitter', 'AdS', 'Anti-de Sitter', 'black hole', 
            'Schwarzschild'.
            Default: `spacetimes.FlatSpacetime` of the determined dimension.
        shape: str, CoordinateShape
            (The name of) a coordinate shape that describes the embedding 
            region of the events. 
            Default: `DefaultShape()` of the spacetime object.
        coordinates: np.ndarray
            List of coordinates, a row of coordinates for event to be created.
        dim: int
            Dimension for the default spacetime.
            Default: 2
        Note: A `ValueError` is raised if the dimensions of any of the four 
        optional arguments are not compatible.
        '''
        # initialise base class (Causet):
        super().__init__()
        # initialise dimension:
        if dim <= 0:
            if (spacetime is not None) and isinstance(spacetime, Spacetime):
                dim = spacetime.Dim
            elif (shape is not None) and isinstance(shape, CoordinateShape):
                dim = shape.Dim
            elif coordinates is not None:
                dim = len(coordinates[0])
            else:
                dim = 2  # default
        # initialise spacetime:
        M: Spacetime
        if spacetime is None:
            M = spacetimes.FlatSpacetime(dim)
        elif isinstance(spacetime, str):
            if spacetime in {'flat', 'Minkowski'}:
                M = spacetimes.FlatSpacetime(dim)
            elif spacetime in {'dS', 'de Sitter'}:
                M = spacetimes.deSitterSpacetime(dim)
            elif spacetime in {'AdS', 'Anti-de Sitter'}:
                M = spacetimes.AntideSitterSpacetime(dim)
            elif spacetime in {'black hole', 'Schwarzschild'}:
                M = spacetimes.BlackHoleSpacetime(dim)
            else:
                raise ValueError(
                    'The spacetime name "%s" is not supported.' % spacetime)
        else:
            M = spacetime
        if M.Dim != dim:
            self.__raiseDimValueError__('spacetime')
        self._spacetime = M
        defaultShape = M.DefaultShape()
        # initialise shape:
        _shape: CoordinateShape
        if shape is None:
            _shape = defaultShape
        elif isinstance(shape, str):
            _shape = CoordinateShape(dim, shape)
        else:
            _shape = shape
        if _shape.Dim != defaultShape.Dim:
            self.__raiseDimValueError__('shape')
        self._shape = _shape
        # create new events:
        if coordinates is not None:
            # add labelled events with coordinates:
            if isinstance(coordinates, np.ndarray) and \
                    ((coordinates.ndim != 2) or
                     (coordinates.shape[1] != defaultShape.Dim)):
                self.__raiseDimValueError__('coordinates')
            self.create(coordinates)

    @staticmethod
    def _Permutation_Coords(P: List[int], radius: float) -> np.ndarray:
        '''
        Returns a matrix of (t, x) coordinates with `len(P)` rows, a pair of 
        coordinates for each element in the permutation integer list (integers 
        from 1 to `len(P)`).
        '''
        count: int = len(P)
        coords: np.ndarray = np.empty((count, 2))
        if count > 0:
            cellscale: float = radius / float(count)
            for i, p in enumerate(P):
                j: int = p - 1
                crd_u: float = (p - 0.5) * cellscale
                crd_v: float = (i + 0.5) * cellscale
                coords[j, 0] = crd_u + crd_v - radius
                coords[j, 1] = crd_u - crd_v
        return coords

    @staticmethod
    def FromPermutation(P: List[int], labelFormat: Optional[str] = None,
                        radius: float=1.0) -> 'EmbeddedCauset':
        '''
        Generates a causal set from the permutation P of integers from 1 to 
        `len(P)` - that can be embedded in an Alexandrov subset of Minkowski 
        spacetime.
        '''
        C: EmbeddedCauset = EmbeddedCauset(
            shape=CoordinateShape(2, 'diamond', radius=radius))
        C.create(EmbeddedCauset._Permutation_Coords(P, radius), labelFormat)
        return C


    ####################################################################
    #___________________________________________________________________
    # EMBEDDING CAUSET FASTER LABELLING
    #___________________________________________________________________
    ####################################################################
    def nlabel (self, reverse = False) -> List[CausetEvent]:
        """
        Turn Causet into list of events with a natural labelling.
        """
        Clist = list(self)
        a = 1 if not reverse else -1      
        Clist.sort(key = lambda e : a*e.Coordinates[0])
        return Clist
    
    def nlist (self, reverse = False) -> List[CausetEvent]:
        """
        Turn Causet into list of events with a natural labelling.
        """
        Clist = list(self)
        a = 1 if not reverse else -1      
        Clist.sort(key = lambda e : a*e.Coordinates[0])
        return Clist


    ####################################################################
    #___________________________________________________________________
    # EMBEDDING PROPERTIES
    #___________________________________________________________________
    ####################################################################
    @property
    def Dim(self) -> int:
        '''
        Returns the coordinate dimension of the embedding region.
        '''
        return self.Shape.Dim

    @property
    def Density(self) -> float:
        '''
        Returns the density of events as ratio of the set cardinality to the 
        embedding shape's volume.
        '''
        return float(self.Card) / self.Shape.Volume

    @property
    def LengthScale(self) -> float:
        '''
        Returns the fundamental length scale as inverse d-root of 
        `self.Density` if `card > 0`, else 0.0. 
        '''
        return 0.0 if self.Card == 0 else self.Density**(1.0 / self.Shape.Dim)

    @property
    def Shape(self) -> CoordinateShape:
        '''
        Returns the CoordinateShape object of the embedding region.
        '''
        return self._shape

    @property
    def Spacetime(self) -> Spacetime:
        '''
        Returns the Spacetime object of the embedding.
        '''
        return self._spacetime


    ####################################################################
    #___________________________________________________________________
    # KINEMATICS
    #___________________________________________________________________
    ####################################################################
    def ptime (self, a: CausetEvent, b: CausetEvent):
        """
        Return (oriented) proper time between element a and b. 
        Note: returns complex value is spacelike interval.
        Not sure if correct in non-flat spacetime.
        """
        avec = a.Coordinates
        bvec = b.Coordinates
        return self._spacetime.ds(avec, bvec)

    def ptime_est(self, a: CausetEvent, b: CausetEvent) -> int:
        '''
        Compute the estimate of proper time between 2 elements a and b.
        '''
        return "This is current a placeholder"
    
    def timegeodesics(self, a: CausetEvent, b: CausetEvent, 
                    mode = "flat"):
        return "This is current a placeholder"

    def MMdim_est(self, d0 = 2, Nsamples = 20,
                ptime_constr = None, 
                size_min = 20,
                optimizer = fsolve,
                opt_sol_index = 0,
                opt_flag_index = 2,
                **optkwargs) -> float:
        '''
        Use Myrheim-Meyers dimensional estimator to compute the 
        fractal dimension (not necesseraly int).

        Parameters:
        --------------------------
        - d0: float\n
            Initial guess for dimension.\n
            Default is 2.

        - Nsamples: int\n
            Number of times to iterate procedure to average on.\
            Default is 20.

        - ptime_constr: callable(ptime) -> Bool\n
            A callable setting a constraint on the proper time between
            the two elements and returning a Boolean.\n
            Default: None.

        - size_min: int\n
            Minimum size of Alexandrov Sets on which apply estimators.\n
            Default is 20 elements.
        
        - optimizer: callable(funct, guess, args, ....) -> list\n
            A function which takes as args: 
            - funct whose roots are to be found
            - initial guess
            - arg to be given for funct
            and returns a array-like obj with:
            - numerical result
            - a flag saying whether it was successful.
            - maybe other stuff\n
            Default is sp.optimizer.fsolve.
        
        - opt_sol_index: int\n
            The optimizer returns a list of stuff. opt_sol_index
            is the index at which the numerical solution is found.\n
            Default 0 (which is the case for most optimizers).
        
        - opt_flag_index: int\n
            The optimizer returns a list of stuff. opt_flag_index
            is the index at which the flag specifying whether it was
            successful or not is found.\n
            Default 2.
        
        - **optkwargs: \n
            kwargs to be passed to optimizer.
        
        Return
        --------------------------
        - dimension estimate: float
        - dimension std: float
        '''

        if not isinstance(self._spacetime, spacetimes.FlatSpacetime):
            print("NOTE for MM Estimator: this is currently working\
                as if the spacetime was flat.")
           
        def MM_drelation(d):
            a = sp.special.gamma(d+1)
            b = sp.special.gamma(d/2)
            c = 4 * sp.special.gamma(3*d/2)
            return a*b/c
        
        def MM_to_solve(d, ord_fr):
            return MM_drelation(d) - ord_fr/2

        N = len(self)
        destimates = []
        isample = 0
        counts = [0,0] #times it worked vs times it did not
        while isample < Nsamples:
            if counts[1]-Nsamples>=1000 and counts[0]==0:
                print("The algorithm never found a suitable Alexandrov\n\
                Interval, whereas it found 1000 unsuitable: probably\n\
                you picked a too small causet. Returning [NaN, NaN]")
                return np.nan, np.nan
            # pick two random elements
            a, b = random.sample(self._events, k = 2)
            # if linked, ptime respects the constraint, 
            # and the size of the interval is big enough
            if CausetEvent.isCausalTo(a, b):
                if ptime_constr is None:
                    if self.IntervalCard(a,b) >= size_min:
                        fr_i = self.ord_fr(self.Interval(a,b, 
                                                disjoin=True))
                        if fr_i == 1:
                            destimates.append(1)
                        else:
                            sol_i = optimizer(MM_to_solve, d0, fr_i,
                                                **optkwargs)
                            if not (opt_flag_index is None):
                                if sol_i[opt_flag_index] == 1: 
                                    d_i= sol_i[opt_sol_index]
                                    destimates.append(d_i)
                                    isample += 1
                                else:
                                    continue
                            else:
                                d_i= sol_i[opt_sol_index]
                                destimates.append(d_i)
                                isample += 1
                elif self.IntervalCard(a,b) >= size_min:
                    counts[0] += 1
                    if ptime_constr(self.ptime(a, b)):
                #Note: switched order between intervalcard and ptime
                # as in generic spacetime ptime might take more.
                        fr_i = self.ord_fr(self.Interval(a,b,
                                        disjoin = True))
                        if fr_i == 1:
                            destimates.append(1)
                        else:
                            sol_i = optimizer(MM_to_solve, d0, fr_i,
                                                **optkwargs)
                            if not (opt_flag_index is None):
                                if sol_i[opt_flag_index] == 1: 
                                    d_i= sol_i[opt_sol_index]
                                    destimates.append(d_i)
                                    isample += 1
                                else:
                                    continue
                            else:
                                d_i= sol_i[opt_sol_index]
                                destimates.append(d_i)
                                isample += 1
                else:
                    counts[1] += 1
                    continue
            else:
                continue

        return np.mean(destimates), np.std(destimates)


    ####################################################################
    #___________________________________________________________________
    # OTHERS
    #___________________________________________________________________
    ####################################################################
    def create(self, coordinates: Union[List[List[float]],
                                        List[np.ndarray],
                                        np.ndarray],
               labelFormat: Optional[str] = None, relate: bool = True) -> \
            Set[CausetEvent]:
        '''
        Creates new events with the specified coordinates, adds them to 
        this instance and returns the new set of events.
        The argument 'coordinates' has to be List[List[float]], 
        List[np.ndarray] or np.ndarray (matrix with a coordinate row for 
        each event).
        '''
        n: int = self.Card + 1
        eventSet: Set[CausetEvent] = {
            CausetEvent(label=(n + i) if labelFormat is None else
                        labelFormat.format(n + i), coordinates=c)
            for i, c in enumerate(coordinates)}
        self._events.update(eventSet)
        if relate:
            self.relate()
        return eventSet

    def relate(self, link: bool = True) -> None:
        '''
        Resets the causal relations between all events based on their 
        embedding in the given spacetime manifold.
        '''
        _iscausal: Callable[[np.ndarray, np.ndarray],
                            Tuple[bool, bool]] = self._spacetime.Causality()
        for e in self._events:
            e._prec = set()
            e._succ = set()
        eventList: List[CausetEvent] = list(self._events)
        eventList_len: int = len(eventList)
        for i, a in enumerate(eventList):
            for j in range(i + 1, eventList_len):
                b = eventList[j]
                isAB, isBA = _iscausal(np.array(a.Coordinates),
                                       np.array(b.Coordinates))
                if isAB:  # A in the past of B:
                    a._succ.add(b)
                    b._prec.add(a)
                if isBA:  # A in the future of B:
                    b._succ.add(a)
                    a._prec.add(b)
        if link:
            self.link()
        else:
            self.unlink()

    def relabel(self, dim: int = 0, descending: bool = False) -> None:
        '''
        Resets the labels of all events to ascending (default) or descending 
        integers (converted to str) corresponding to the coordinate component 
        in dimension 'dim'.
        '''
        eventList = list(self._events)
        sorted_idx = np.argsort(np.array(
            [e.Coordinates[dim] for e in eventList]))
        if descending:
            sorted_idx = np.flip(sorted_idx)
        for i, idx in enumerate(sorted_idx):
            eventList[idx].Label = i + 1
