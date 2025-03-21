import numpy as np
from scipy import constants as spc
import optical_depth

def EUV_flux(z_0: float, solar_zenith_angle: float, absorption_crosssecs: np.ndarray, 
                  number_densities: np.ndarray, heights: np.ndarray,  irradiance: float, wavelength: float, 
                  planet_radius: float = 6378e3) -> float:
    
    '''
    Function that calculates the solar EUV photon flux at a given altitude using the optical depth and
    irradiance data

    Parameters:
        z_0 (float): Altitude [m] at which to compute the flux.
        solar_zenith_angle (float): Solar zenith angle [degrees].
        absorption_crosssecs (np.ndarray): Absorption cross-sections [m^2] for relevant species.
        number_densities (np.ndarray): Number densities [m^-3] of atmospheric species at all altitudes.
        heights (np.ndarray): Altitudes [m] corresponding to number densities.
        irradiance (float): Unattenuated spectral irradiance [W/m^2/m] at the top of the atmosphere.
        wavelength (float): Wavelength [m] corresponding to the irradiance value.
        planet_radius (float, optional): Planetary radius [m]. Default is Earth's mean radius (6378e3 m).

    Returns:
        float: Photon flux [photons/m^2/s/nm]
    '''

    #FInding the optical depth for the given parameters
    tau = optical_depth.optical_depth(z_0, solar_zenith_angle, absorption_crosssecs,
                                      number_densities, heights, planet_radius)
    
    #Defining constants with scipy
    c = spc.speed_of_light
    h = spc.h

    #Converting to photon flux instead of energy
    inf_flux = irradiance * (wavelength / (h * c))

    #Calculating the flux
    flux = inf_flux * np.exp(- tau)

    return flux