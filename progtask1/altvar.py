import numpy as np
import scipy.constants as spc

def altitude_var_g(heights: np.ndarray) -> np.ndarray:
    '''Calculates the gravitational acceleration for given heights above ground.
    Takes in an array of heights above ground in meters and returns array
    of gravitational acceleration in m/s^2'''
    
    #Defining the constants, Earth radius in meters, gravitational acceleration at ground in m/s^2
    earth_radius = 6378000 #m
    grav_acceleration_at_ground = spc.g #m/s^2

    #Using the given formula to calculate gravitational acceleration at each height
    grav_acc = grav_acceleration_at_ground * (earth_radius**2 / (earth_radius + heights)**2)

    return grav_acc

#testing...
if __name__ == '__main__':

    heigths_in_m = np.arange(0, 1000000, 1000) #array containing heights with ubit meters
    print(altitude_var_g(heights=heigths_in_m))

