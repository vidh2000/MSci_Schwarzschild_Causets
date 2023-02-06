import numpy as np


# Areas in l^2 units
#areas_ll = 4*np.pi*(2*mass)**2/(scale**2)

areas_ll = np.arange(1000,20001,1000)
rhos = np.array([5000])
int8_t = False

for Rho in rhos:
    print("")
    
    # Set scale (discreteness length)
    scale = Rho **(-1/4)
    print(f"Rho={Rho:.0f}, Scale = {scale:.3f}")
   
    masses = np.sqrt(areas_ll/(16*np.pi*np.sqrt(Rho)))
    print("Masses:")
    print(masses, "\n")

    for mass in masses:

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
        print(f"Rho={Rho}, M={mass:.2f}, Area={4*np.pi*(2*mass)**2/(scale**2):.0f} l^2, N={N:.0f} \
        RAM={mem:.0f}GB") #add time later
