import numpy as np

mass = 1
Rho = np.array([400, 500, 600, 700, 800, 900, 1000])
scale = Rho **(-1/4)
R = 2*mass+3*scale
r = 2*mass-3*scale
T = 4*scale
h = r/R
N = Rho * (4*3.1415/3) * (R*R*R-r*r*r) * T
print("\n",N)

mass = 2
Rho = np.array([400, 500, 600, 700, 800, 900, 1000])
scale = Rho **(-1/4)
R = 2*mass+3*scale
r = 2*mass-3*scale
T = 4*scale
h = r/R
N = Rho * (4*3.1415/3) * (R*R*R-r*r*r) * T
print("\n",N)

mass = 3
Rho = np.array([400, 500, 600, 700, 800, 900, 1000, 1200])
scale = Rho **(-1/4)
R = 2*mass+3*scale
r = 2*mass-3*scale
T = 4*scale
h = r/R
N = Rho * (4*3.1415/3) * (R*R*R-r*r*r) * T
print("\n",N)