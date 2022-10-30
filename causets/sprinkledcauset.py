#!/usr/bin/env python
'''
Created on 22 Jul 2020

@author: Christoph Minz
@license: BSD 3-Clause
'''
from __future__ import annotations
from typing import Set, List, Iterable, Union, Optional
import numpy as np
import random as rand
from numpy.random import default_rng
from causets.causetevent import CausetEvent  # @UnresolvedImport
from causets.embeddedcauset import EmbeddedCauset  # @UnresolvedImport
from causets.shapes import CoordinateShape  # @UnresolvedImport
from causets.spacetimes import Spacetime  # @UnresolvedImport
from copy import deepcopy

class SprinkledCauset(EmbeddedCauset):
    '''
    Handles a causal set that is embedded in a subset of a manifold.
    '''

    _intensity: float

    def __init__(self, card: int = 0, intensity: float = 0.0, dim: int = -1,
                 spacetime: Optional[Spacetime] = None,
                 shape: Optional[Union[str, CoordinateShape]] = None) -> None:
        '''
        Generates a sprinkled causal set by sprinkling in a spacetime subset. 
        The arguments `dim`, `shape` and `spacetime` are handled by the super 
        class `EmbeddedCauset` before events are sprinkled.

        Parameters
        ----------
        card: int
            Number of sprinkled events.
        intensity: float
            Sprinkling intensity parameter, the expected number of sprinkled 
            events.
        '''
        # initialise base class (EmbeddedCauset):
        super().__init__(spacetime=spacetime, shape=shape, dim=dim)
        # sprinkle:
        self._intensity = 0.0
        if card > 0:
            self.sprinkle(card)
        else:
            self.intensify(intensity)
         
    @property
    def Intensity(self) -> float:
        '''
        Returns the sprinkling intensity, which is the expected number of 
        sprinkled events. The exact number of sprinkled events is given by the 
        property 'Card'.
        '''
        return self._intensity

    @property
    def Density(self) -> float:  # overwrites superclass
        return self._intensity / self.Shape.Volume

    @property
    def LengthScale(self) -> float:  # overwrites superclass
        return (self.Shape.Volume / self._intensity)**(1.0 / self.Dim)

    def _sprinkle_coords(self, count: int, shape: CoordinateShape, rng) -> \
            np.ndarray:
        if count < 0:
            raise ValueError(
                'The sprinkle cardinality has to be a non-negative integer.')
        coords: np.ndarray = np.empty((count, self.Dim), dtype=np.float64)
        if shape.Name in ('cube', 'cuboid'):
            # Create rectangle based sprinkle:
            low: np.ndarray
            high: np.ndarray
            if shape.Name == 'cuboid':
                low = shape.Center - shape.Parameter('edges') / 2
                high = shape.Center + shape.Parameter('edges') / 2
            else:
                low = shape.Center - shape.Parameter('edge') / 2
                high = shape.Center + shape.Parameter('edge') / 2
            for i in range(count):
                coords[i, :] = rng.uniform(low, high)
        elif shape.Name in ('ball', 'cylinder', 'diamond', "bicone"):
            # Create circle based sprinkle:
            isCylindrical: bool = 'cylinder' in shape.Name
            isDiamond: bool = 'diamond' in shape.Name or "bicone" in shape.Name

            d: int = self.Dim
            radius: float = shape.Parameter('radius')
            if (d == 2) and isDiamond:
                
                # pick `count` random coordinate tuples uniformly:
                uv: np.ndarray = np.random.uniform(low=-1.0, high=1.0,
                                                    size=(count, 2))
                coords[:, 0] = uv[:, 0] + uv[:, 1]
                coords[:, 1] = uv[:, 0] - uv[:, 1]
                coords *= radius / 2
            else:
                # n_rad is, at least de facto, number of "radial" coordinates
                rad_start: int = 0 if shape.Name == 'ball' else 1
                n_rad: int = d - rad_start
                if isCylindrical:
                    # set time coordinate:
                    time_low: float = shape.Center[0] - \
                        shape.Parameter('duration') / 2
                    time_high: float = shape.Center[0] + \
                        shape.Parameter('duration') / 2
                    coords[:, 0] = rng.uniform(time_low, time_high,
                                               size=(count,))
                # pick `count` random coordinate tuples uniformly:
                r_low: float = shape.Parameter('hollow')**n_rad
                for i in range(count):
                    # EXPLAINED at
                    # https://math.stackexchange.com/a/87238
                    # DERIVED from
                    # Muller, M. E. "A Note on a Method for Generating Points 
                    # Uniformly on N-Dimensional Spheres." Comm. Assoc. Comput. 
                    # Mach. 2, 19-20, Apr. 1959.
                    # and
                    # Marsaglia, G. "Choosing a Point from the Surface of 
                    # a Sphere." Ann. Math. Stat. 43, 645-646, 1972.
                    # note: all for ball, spatials for cylinder/bicone
                    # 1) Get coordinates using normal distribution
                    coord: np.ndarray = rng.standard_normal(size=(n_rad,))
                    # 2) get the radius associated to these
                    r: float = np.sqrt(sum(np.square(coord)))
                    # 3) normalise: to uniformly distribute on a sphere of r=1
                    coord /= r
                    # 4) Get scaling deviate (note lower boundary if hollow)
                    r_scaling: float = rng.uniform(low=r_low)**(1.0 / n_rad)
                    # 5) Scale Coordinates
                    coord *= radius * r_scaling
                    if isDiamond:
                        # need to set time coordinate properly
                        # 1) take d-root of uniform deviate
                        droot_x: float = rng.uniform()**(1.0 / d)
                        # 2) pick if it goes in upper or lower cone
                        t_sign: float=np.sign(rng.uniform(low=-1.0, high=1.0))
                        # 3) apply trasnformation method: uniform -> time PDF
                        coords[i, 0] = t_sign * (1 - droot_x) * radius
                        # 4) Adjust scaling of radius based on t
                        coord *= droot_x
                    coords[i, rad_start:] = shape.Center[rad_start:] + coord
        return coords


    def sprinkle(self, count: int, rng=default_rng(),
                 shape: Optional[CoordinateShape] = None) -> Set[CausetEvent]:
        '''
        Creates a fixed number of new events by sprinkling into `shape` (by 
        default the entire embedding region).
        '''
        if count < 0:
            raise ValueError(
                'The sprinkle cardinality has to be a non-negative integer.')
        self._intensity += float(count)
        if shape is None:
            shape = self.Shape
        coords: np.ndarray = self._sprinkle_coords(count, shape, rng)
        return super().create(coords)

    def intensify(self, intensity: float, rng=default_rng(),
                  shape: Optional[CoordinateShape] = None) -> Set[CausetEvent]:
        '''
        Creates an expected number of new events by sprinkling into `shape` (by 
        default the entire embedding region). The expected number is determined 
        by the Poisson distribution with the given `intensity` parameter.
        '''
        if intensity < 0.0:
            raise ValueError(
                'The intensity parameter has to be a non-negative float.')
        self._intensity += intensity
        count: int = int(rng.poisson(lam=intensity))
        if shape is None:
            shape = self.Shape
        coords = self._sprinkle_coords(count, shape, rng)
        return super().create(coords)

    def create(self, coords: Union[Iterable[List[float]],
                                   Iterable[np.ndarray],
                                   np.ndarray],
               labelFormat: Optional[str] = None, relate: bool = True) -> \
            Set[CausetEvent]:
        card_old: float = float(self.Card)
        eventSet: Set[CausetEvent] = super().create(
            coords, labelFormat, relate)
        self._intensity += (float(self.Card) - card_old)
        return eventSet

    def add(self, eventSet: Iterable, unlink: bool = False) -> None:
        card_old: float = float(self.Card)
        super().add(eventSet, unlink)
        self._intensity += (float(self.Card) - card_old)

    def discard(self, eventSet: Iterable, unlink: bool = False) -> None:
        card_old: float = float(self.Card)
        super().discard(eventSet, unlink)
        self._intensity *= (float(self.Card) / card_old)
