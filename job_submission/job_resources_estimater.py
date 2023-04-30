import numpy as np


# Areas in l^2 units
#areas_ll = 4*np.pi*(2*mass)**2/(scale**2)

dim=4
areas_ll = np.arange(1000,30001,500)
rhos = np.array([5000])
int8_t = False

for Rho in rhos:
    print("")
    
    # Set scale (discreteness length)
    if dim==4:
        scale = Rho **(-1/4)
        #masses = np.sqrt(areas_ll/(16*np.pi*np.sqrt(Rho))) #equiv to line below
        masses = np.sqrt(areas_ll/(16*np.pi)*scale**2)
        print("Masses:")
        print(masses, "\n")
        print([round(mass,2) for mass in masses])

    elif dim==3:
        scale = Rho**(-1/3)
        print(f"Rho={Rho:.0f}, Scale = {scale:.3f}")
        masses = areas_ll*scale/(4*np.pi)
        print("Masses:")
        print(masses, "\n")
        print([round(mass,2) for mass in masses])
   

    for mass in masses:

        # Set sprinkling region
        R = 2*mass+3*scale
        r = 2*mass-3*scale
        T = 4*scale

        # Calculate the cardinality of the sprinkled causet
        if dim==4:
            N = Rho * (4*3.1415/3) * (R*R*R-r*r*r) * T
        elif dim==3:
            N = Rho * 3.1415 * (R*R-r*r) * T    
        # Find required RAM
        if int8_t:
            mem = N**2 / 1073741824   #in GB
        else:
            mem = N**2 * 4 / 1073741824 #in GB


        
        # Print require information
        if dim==4:
            print(f"Rho={Rho}, M={mass:.2f},\
    Area={4*np.pi*(2*mass)**2/(scale**2):.0f} l^2, N={N:.0f} \
        RAM={mem:.0f}GB") #add time later

        elif dim==3:
            print(f"Rho={Rho}, M={mass:.2f},\
    Area={2*np.pi*(2*mass)/(scale):.0f} l, N={N:.0f} \
        RAM={mem:.0f}GB") #add time later


# print("\nGET N FOR RHO=",rhos[0]," and M=1")
# mass = 1
# scale = rhos[0] **(-1/4)
# R = 2*mass+3*scale
# r = 2*mass-3*scale
# T = 4*scale
# N = rhos[0] * (4*3.1415/3) * (R*R*R-r*r*r) * T
# mem = N**2 * 4 / 1073741824 #in GB
# print(f"Area={4*np.pi*(2*mass)**2/(scale**2):.0f} l^2, N={N:.0f} \
#     RAM={mem:.0f}GB")
