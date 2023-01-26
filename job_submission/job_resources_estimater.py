import numpy as np


masses = [1, 1.5, 2, 2.5, 2.6, 3]
rhos = np.array([5000])
int8_t = False

for Rho in rhos:
    print("")
    for mass in masses:
        
        # Set scale (discreteness length)
        scale = Rho **(-1/4)

        # Set sprinkling region
        R = 2*mass+3*scale
        r = 2*mass-3*scale
        T = 4*scale

        # Calculate the cardinality of the sprinkled causet
        N = Rho * (4*3.1415/3) * (R*R*R-r*r*r) * T
        
        # Find required RAM
        if int8_t:
            mem = N**2 / 1073741824   #in GB
        else:
            mem = N**2 * 4 / 1073741824 #in GB

        # Print require information
        print(f"Rho={Rho}, M={mass:.2f}, N={N:.0f} \
        RAM={mem:.0f}GB") #add time later

