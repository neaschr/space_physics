import numpy as np

#We will make functions to find the optical depth as a function of altitude and wavelength

def optical_depth(z_0: float, solar_zenith_angle: float, absorption_crosssecs: np.ndarray, 
                  number_densities: np.ndarray, heights: np.ndarray, planet_radius: float = 6378e3) -> float:
    
    '''
    Compute the optical depth using numerical integration.

    Parameters:
    - z_0 (float): Starting altitude (in meters)
    - solar_zenith_angle (float): Solar zenith angle (in degrees)
    - absorption_crosssecs (np.ndarray): Absorption cross-section for each particle species 
    - number_densities (np.ndarray): Number densities by height for each particle species (in m-3)
    - heights (np.ndarray): Heights corresponding to the number densities (in meters)
    - planet_radius (float): Planet radius in meters (default earth radius: 6378 km)

    Returns:
    - float: Optical depth for one wavelength
    '''
    #Convert solar zenith angle to radians
    solar_zen_angle_rad = np.radians(solar_zenith_angle)

    #Find the index corresponding to z_0
    z0_idx = np.searchsorted(heights, z_0)

    #We only want the heights and number densities above z_0
    number_densities = number_densities[:, z0_idx:]  
    heights = heights[z0_idx:]  
    
    if solar_zenith_angle < 70: 
        
        #Numerically integrating the number densities by height (given in meters), using trapezoidal integration
        integral = np.trapz(number_densities, heights) 

        #Defining the constant outside the sum
        const = (1 / (np.cos(solar_zen_angle_rad)))
        
        #Finding the optical depth
        opt_depth = const * np.sum(absorption_crosssecs * integral, axis=0)

    elif solar_zenith_angle >= 70 and solar_zenith_angle < 90:

        #Finding the correction to the number densities
        factor = (1 - ((planet_radius + z_0) / (planet_radius + heights))**2 * (np.sin(solar_zen_angle_rad))**2)**(-(1/2))
        number_densities2 = number_densities * factor  

        #Numerically integrating the number densities by height (given in meters), using trapezoidal integration
        integral = np.trapz(number_densities2, heights)

        #Finding the optical depth
        opt_depth = np.sum(absorption_crosssecs * integral, axis=0)
        
    else:
        opt_depth = 0
        print('Nea is sad')

    return opt_depth

