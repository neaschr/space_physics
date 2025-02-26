import numpy as np
from scipy import constants as spc
import matplotlib.pyplot as plt
import datareader

#We will make functions to find the optical depth as a function of altitude and wavelength

def optical_depth(z_0: float, solar_zenith_angle: float, absorption_crosssecs: np.ndarray, number_densities: np.ndarray, heights: np.ndarray, planet_radius: float = 6378e3):
    
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
    - np.ndarray: Optical depth for one wavelength
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

if __name__ == '__main__':

    #Importing the data for absorption crossections
    names_of_columns = [
        'wavelength [m]', 
        'absorption cross section N2 [m^2]',
        'absorption cross section O  [m^2]', 
        'absorption cross section O2 [m^2]'
    ]

    absorp_crosssec_data = datareader.data_reader('phot_abs.dat', names_of_columns)

    names_of_columns2 = [
        'Height_km',
        'Oksygen-atoms_cm-3',
        'Nitrogen-molecules_cm-3',
        'Oksygen-molecules_cm-3',
        'Mass_density_g_per_cm3',
        'Temperature_neutral_K'
    ]

    data = datareader.data_reader('MSIS.dat', names_of_columns2)

    #Defining arrays for the number density of each species, number density in n/cm^3
    N_atom_oxy = data['Oksygen-atoms_cm-3'].to_numpy() 
    N_molec_oxy = data['Oksygen-molecules_cm-3'].to_numpy()
    N_molec_nitro = data['Nitrogen-molecules_cm-3'].to_numpy()

    #Converting to per meter cubed
    N_atom_oxy = N_atom_oxy * 1e6
    N_molec_oxy = N_molec_oxy * 1e6
    N_molec_nitro = N_molec_nitro * 1e6
    
    #Defining array for the heights, given in unit meters
    height_in_meters = data['Height_km'].to_numpy() * 1e3

    optical_depth_calculated = np.zeros(len(height_in_meters))

