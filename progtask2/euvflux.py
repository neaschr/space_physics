import numpy as np
import optical_depth

def EUV_flux(z_0: float, solar_zenith_angle: float, absorption_crosssecs: np.ndarray, 
                  number_densities: np.ndarray, heights: np.ndarray, planet_radius: float = 6378e3) -> float:
    
    tau = optical_depth.optical_depth(z_0, solar_zenith_angle, absorption_crosssecs,
                                      number_densities, heights, planet_radius)
    
    inf_flux = 10e10
    flux = inf_flux * np.exp(- tau)

    return flux



