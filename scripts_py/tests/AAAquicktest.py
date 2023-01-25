import numpy as np


masses = [1,2,3]
Rho = np.array([1000,5000,10000])

print("\n Rhos =", Rho)

for mass in masses:

    scale = Rho **(-1/4)
    R = 2*mass+3*scale
    r = 2*mass-3*scale
    T = 4*scale
    h = r/R
    N = Rho * (4*3.1415/3) * (R*R*R-r*r*r) * T
    print(f"M = {mass}")
    print([int(n) for n in N])

