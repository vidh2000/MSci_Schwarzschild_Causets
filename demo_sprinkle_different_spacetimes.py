#!/usr/bin/env python
'''
Created on 9 Oct 2020

@author: Christoph Minz
@license: BSD 3-Clause
'''

from __future__ import annotations
from typing import List, Tuple  # @UnusedImport
from causets.sprinkledcauset import SprinkledCauset  
from causets.spacetimes import *
from causets.shapes import CoordinateShape 
from causets.causetevent import CausetEvent 
import causets.causetplotting as cplt  
import matplotlib.pyplot as plt



#%% Minkowski 1+1 spacetime
figname = "Minkowski 1+1 spacetime"
S: CoordinateShape = CoordinateShape(2, 'cylinder', duration=2.0, hollow=0)
S: CoordinateShape = CoordinateShape(2, 'cuboid', edges=[3,2])
C: SprinkledCauset = SprinkledCauset(intensity=10.0,
                                     spacetime=FlatSpacetime(2),
                                     #spacetime="Minkowski",
                                     shape=S)
e: CausetEvent = C.CentralAntichain().pop()  # pick one event


plt.figure(figname)
cplt.setDefaultColors('UniYork')  # using University of York brand colours
dims: List[int] = [1, 0]  # choose the (order of) plot dimensions
if len(dims) > 2:
    cplt.figure(figsize=(8.0, 8.0))
S.plot(dims)  # plot the embedding shape
# Add causet plots and show result:
cplt.plot(C, dims=dims, events={'alpha': 1},
          links={'alpha': 1, 'linewidth': 0.3}, labels=False)
#cplt.plot(list(e.Cone), dims=dims, spacetime=C.Spacetime,
#          events={'markerfacecolor': 'cs:darkblue'},
#          links={'alpha': 0.6, 'linewidth': 1.5}, labels=False)
# end times for past and future light-cone:
#timeslices: Tuple[float, float] = S.Limits(0)
#cplt.plot(e, dims=dims, spacetime=C.Spacetime,
#          events={'markerfacecolor': 'cs:red'},
#          pastcones={'alpha': 1.0}, futurecones={'alpha': 1.0},
#          time=timeslices)
ax: cplt.Axes = cplt.gca()
ax.set_xlabel('space' if dims[0] > 0 else 'time')
ax.set_ylabel('space' if dims[1] > 0 else 'time')
if len(dims) > 2:
    ax.set_zlabel('space' if dims[2] > 0 else 'time')
    ax.grid(False)
cplt.show()



C = list(C)
print(C)

#%% Minkowski 2+1 spacetime


S: CoordinateShape = CoordinateShape(3, 'cylinder', duration=2.0, hollow=0)
#S: CoordinateShape = CoordinateShape(3, 'cuboid', edges=[3,2,2])
C: SprinkledCauset = SprinkledCauset(intensity=100.0,
                                     spacetime=FlatSpacetime(3),
                                     #spacetime="Minkowski",
                                     shape=S)
e: CausetEvent = C.CentralAntichain().pop()  # pick one event


plt.figure("Minkowski 2+1 spacetime")
cplt.setDefaultColors('UniYork')  # using University of York brand colours
dims: List[int] = [1,2,0]  # choose the (order of) plot dimensions
if len(dims) > 2:
    cplt.figure(figsize=(8.0, 8.0))
S.plot(dims)  # plot the embedding shape
# Add causet plots and show result:
cplt.plot(C, dims=dims, events={'alpha': 0.05},
          links={'alpha': 1, 'linewidth': 0.3}, labels=False)
#cplt.plot(list(e.Cone), dims=dims, spacetime=C.Spacetime,
#          events={'markerfacecolor': 'cs:darkblue'},
#          links={'alpha': 0.6, 'linewidth': 1.5}, labels=False)
# end times for past and future light-cone:
#timeslices: Tuple[float, float] = S.Limits(0)
#cplt.plot(e, dims=dims, spacetime=C.Spacetime,
#          events={'markerfacecolor': 'cs:red'},
#          pastcones={'alpha': 1.0}, futurecones={'alpha': 1.0},
#          time=timeslices)
ax: cplt.Axes = cplt.gca()
ax.set_xlabel('space' if dims[0] > 0 else 'time')
ax.set_ylabel('space' if dims[1] > 0 else 'time')
if len(dims) > 2:
    ax.set_zlabel('space' if dims[2] > 0 else 'time')
    ax.grid(False)
cplt.show()

#%% Black-hole spacetime


S: CoordinateShape = CoordinateShape(2, 'ball', radius=3.0)
S: CoordinateShape = CoordinateShape(2, 'cuboid', edges=[3,3])
C: SprinkledCauset = SprinkledCauset(intensity=500.0,
                                     spacetime=BlackHoleSpacetime(2),
                                     #spacetime="Minkowski",
                                     shape=S)
e: CausetEvent = C.CentralAntichain().pop()  # pick one event


plt.figure("Minkowski 2+1 spacetime")
cplt.setDefaultColors('UniYork')  # using University of York brand colours
dims: List[int] = [1,0]  # choose the (order of) plot dimensions
if len(dims) > 2:
    cplt.figure(figsize=(8.0, 8.0))
S.plot(dims)  # plot the embedding shape
# Add causet plots and show result:
cplt.plot(C, dims=dims, events={'alpha': 0.05},
          links={'alpha': 1, 'linewidth': 0.3}, labels=False)
#cplt.plot(list(e.Cone), dims=dims, spacetime=C.Spacetime,
#          events={'markerfacecolor': 'cs:darkblue'},
#          links={'alpha': 0.6, 'linewidth': 1.5}, labels=False)
# end times for past and future light-cone:
#timeslices: Tuple[float, float] = S.Limits(0)
#cplt.plot(e, dims=dims, spacetime=C.Spacetime,
#          events={'markerfacecolor': 'cs:red'},
#          pastcones={'alpha': 1.0}, futurecones={'alpha': 1.0},
#          time=timeslices)
ax: cplt.Axes = cplt.gca()
ax.set_xlabel('space' if dims[0] > 0 else 'time')
ax.set_ylabel('space' if dims[1] > 0 else 'time')
if len(dims) > 2:
    ax.set_zlabel('space' if dims[2] > 0 else 'time')
    ax.grid(False)
cplt.show()
