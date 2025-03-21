import numpy as np

def irradiance_at_wavelength(wavelength: float, wavelengths: np.ndarray, irradiances: np.ndarray) -> float:
    
    '''
    Function that interpolates irradiance at a given wavelength

    This function performs a linear interpolation to estimate the spectral irradiance
    at a given wavelength using arrays of known wavelengths and their corresponding
    irradiance

    Parameters:
        wavelength (float): The target wavelength [m] to find the irradiance for
        wavelengths (np.ndarray): Array of wavelengths [m]
        irradiances (np.ndarray): Array ofirradiance values [W/m^2/m] corresponding to the wavelengths

    Returns:
        float: Interpolated irradiance [W/m^2/m] at the given wavelength
    '''
    #interpolating to find irradiance
    irradiance =  np.interp(wavelength, wavelengths, irradiances)
    
    return irradiance