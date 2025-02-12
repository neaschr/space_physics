import numpy as np
import scipy.constants as spc

def scale_height(temperature: np.ndarray, particle_mass: np.ndarray, grav_acc=spc.g) -> np.ndarray:
    '''Function that takes in gravitational acceleration, and arrays of particle temperature and particle mass.
    Calculates scale height at each interval'''

    #Defining the constants used, k is the boltzmann constant
    k = spc.k #J/K
    
    #Finding scale height
    scale_h = k * temperature / (particle_mass * grav_acc)
    scale_h = scale_h/1000 #Converting to km
    
    return scale_h

#For 3 divide n for one species by dndz for that species, use numoy gradient