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


    def create_EmbeddedCauset_from_file(self, file):
        """
        Providing coordinates, shape and Spacetime and causal matrix from
        a text/csv file for naturally labeled events it creates the
        EmbeddedCauset object.
        The file either contains cmatrix or pasts+futures.
        This is declared in the first line of the file as the
        "Storage option".

        Parameters inside:
            - storage_option
                        string: Either "sets" or "cmatrix"
            - Size      Number of events sprinkled.
            - Dimension Dimension of the spacetime in which we sprinkled.
            - spacetime string - specifies the spacetime from which spacetime
                        interval is obtained.
                        e.g. "flat", "Minkowski", "black hole", "Schwarzschild"
            - coords    Coordinates vector for N events in dim-D space (N,dim)
            THEN depending on the storage_option parameter:
            - cmatrix   NxN matrix containing of 1 (-1) or 0
                        ->upper triangular matrix required input.
            - past and future sets
                        vector of vectors of size (N,setsize_n)

        """
        storage_option = str(np.genfromtxt(file, delimiter=',', dtype='unicode',
                            usecols=[1])[0])
        self.size = int(np.genfromtxt(file, delimiter=',', dtype='unicode',
                            usecols=[1])[1])
        self.dim = int(np.genfromtxt(file, delimiter=',', dtype='unicode',
                            usecols=[1])[2])
        self.shape_name = str(np.genfromtxt(file, delimiter=',', dtype='unicode',
                            usecols=[1])[3])
        self.spacetime_name = str(np.genfromtxt(file, delimiter=',', dtype='unicode',
                            usecols=[1])[4])


        lines = [[v for v in line.split(",")] for line in open(file)]
        if storage_option == "cmatrix":
            self.cmatrix = np.genfromtxt(file, delimiter=',', dtype='unicode',
                            usecols=[6, 6+self.size-1])
            raise ValueError("Please choose 'sets' storage_option.\
                            'cmatrix' is yet to be implemented properly...")
        elif storage_option == "sets":
            self.pasts = [[int(v) for v in line if v != "\n"] for line in
                                lines[6:6+self.size]]
            self.futures = [[int(v) for v in line if v != "\n"] for line in
                                lines[6+self.size+1:6+2*self.size+1]]
            self.past_links = [[int(v) for v in line if v != "\n"] for line in
                                lines[6+2*self.size+2:6+3*self.size+2]]
            self.fut_links = [[int(v) for v in line if v != "\n"] for line in
                                lines[6+3*self.size+3:6+4*self.size+3]]
        else: 
            raise ValueError("Stored data is not in the format of 'sets' or 'cmatrix'.")
        
        self.coords = [[float(v) for v in line if v != "\n"] for line in
                                            lines[-self.size:]]
    
        # Add relations (futs,pasts,links...) to Events in EventSet
        if storage_option == "sets":
            eventSet: List[CausetEvent] = []
            for i, coord in enumerate(self.coords):
                e = CausetEvent(label=(1 + i), coordinates=coord)
                eventSet.append(e)
            self._events = eventSet
            for (e,past,fut,plink,flink) in zip(self._events, self.pasts,
                            self.futures,self.past_links,self.fut_links):
                e._prec = set(self._events[k] for k in past)
                e._succ = set(self._events[k] for k in fut)
                e._lprec = set(self._events[k] for k in plink)
                e._lsucc = set(self._events[k] for k in flink)
            self._events = set(self._events)
        elif storage_option == "cmatrix":
            raise ValueError("Please choose 'sets' storage_option.\
                            'cmatrix' is yet to be implemented properly...")
      


    
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
    
    def get_coords(self):
        """
        Returns matrix of coordinates, with ith entry being 
        dim-dimensional vector of coordinates of ith event
        """
        return [ei._coordinates for ei in self._events] 


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
        if mode =="flat":
            pass
        return "This is current a placeholder"

    def MMdim_est(self, method = "random", d0 = 2, 
                Nsamples = 20, size_min = 10, size_max = np.Inf,
                ptime_constr = None, 
                optimizer = fsolve,
                opt_sol_index = 0,
                opt_flag_index = 2,
                **optkwargs) -> float:
        '''
        Use Myrheim-Meyers dimensional estimator to compute the 
        fractal dimension (not necesseraly int).

        Parameters:
        --------------------------
        - method: str \n
            - 'random': randomly sample.
            - 'big': take all events with no past, all with no future
            and apply estimator ro their combinations.

        - d0: float \n
            Initial guess for dimension.\n
            Default is 2.

        - Nsamples: int \n
            Times to iterate procedure to then average on if method "random".\
            Default is 20.

        - size_min: int\n
            Minimum size of Alexandrov Sets on which you apply estimators.\n
            Default is 20 elements.
        
        - size_max: int\n
            Maximum size of Alexandrov Sets on which you apply estimators.\n
            Default and highly recommended is np.Inf.

        - ptime_constr: callable(ptime) -> Bool
            A callable setting a constraint on the proper time between
            the two elements and returning a Boolean.\n
            Default: None.

        
        Optimizer-Related Parameters
        -------------------------
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

        if method == "random":
            isample = 0
            fails = 0
            successes = 0

            while isample < Nsamples:
                
                if fails>=1e3 and successes==0:
                    print("The algorithm never found a suitable\
                    Alexandrov Interval, whereas it found 1000 unsuitable.\
                    You picked a too small causet. Returning [NaN, NaN]\n")
                    return np.nan, np.nan
                
                # pick two random elements
                e1, e2 = random.sample(self._events, k = 2)
                if e1 == e2:
                    fails += 1
                    continue
                elif e1 < e2:
                    a = e1
                    b = e2
                elif e1 > e2:
                    a = e2
                    b = e1
                else:
                    fails += 1
                    continue
    
                n = self.IntervalCard(a,b)
                if n >= size_min and n <= size_max:
                    
                    if ptime_constr is None:
                        successes += 1
                        fr_i = self.ord_fr_ab(a,b, den = 'choose')
                        if fr_i == 1:
                            destimates.append(1)
                            isample += 1
                        else:
                            fr_i *= (n-1)/n #correction for MMestimator
                            sol_i = optimizer(MM_to_solve, d0, fr_i,**optkwargs)
                            if not (opt_flag_index is None):
                                if sol_i[opt_flag_index] == 1: 
                                    d_i= sol_i[opt_sol_index][0]
                                    destimates.append(d_i)
                                    isample += 1
                                else:
                                    continue
                            else:
                                d_i = sol_i[opt_sol_index][0]
                                destimates.append(d_i)
                                isample += 1
                
                    elif ptime_constr(self.ptime(a, b)):
                        sucesses += 1
                        fr_i = self.ord_fr_ab(a,b, den = 'choose')
                        if fr_i == 1:
                            destimates.append(1)
                            isample += 1
                        else:
                            fr_i *= (n-1)/n #correction for MMestimator
                            sol_i = optimizer(MM_to_solve, d0, fr_i, **optkwargs)
                            if not (opt_flag_index is None):
                                if sol_i[opt_flag_index] == 1: 
                                    d_i= sol_i[opt_sol_index][0]
                                    destimates.append(d_i)
                                    isample += 1
                                else:
                                    continue
                            else:
                                d_i= sol_i[opt_sol_index][0]
                                destimates.append(d_i)
                                isample += 1
                else:
                    fails += 1
                    continue
        
        elif method == "big":
            As = []
            Bs = []
            for e in self._events:
                if e.PastCard == 0:
                    As.append(e)
                elif e.FutureCard == 0:
                    Bs.append(e)
            
            for i, a in enumerate(As):
                for j, b in enumerate (Bs):
                    
                    n = self.IntervalCard(a,b)
                    if n >= size_min and n <= size_max:
                        
                        if ptime_constr is None:
                            fr_i = self.ord_fr_ab(a,b, den = 'choose')
                            if fr_i == 1:
                                destimates.append(1)
                            else:
                                fr_i *= (n-1)/n #correction for MMestimator
                                sol_i = optimizer(MM_to_solve, d0, fr_i,
                                                  **optkwargs)
                                if not (opt_flag_index is None):
                                    if sol_i[opt_flag_index] == 1: 
                                        d_i= sol_i[opt_sol_index][0]
                                        destimates.append(d_i)
                                    else:
                                        continue
                                else:
                                    d_i = sol_i[opt_sol_index][0]
                                    destimates.append(d_i)
                    
                        elif ptime_constr(self.ptime(a, b)):
                            fr_i = self.ord_fr_ab(a,b, den = 'choose')
                            if fr_i == 1:
                                destimates.append(1)
                            else:
                                fr_i *= (n-1)/n #correction for MMestimator
                                sol_i = optimizer(MM_to_solve, d0, fr_i,
                                                 **optkwargs)
                                if not (opt_flag_index is None):
                                    if sol_i[opt_flag_index] == 1: 
                                        d_i= sol_i[opt_sol_index][0]
                                        destimates.append(d_i)
                                    else:
                                        continue
                                else:
                                    d_i= sol_i[opt_sol_index][0]
                                    destimates.append(d_i)
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


