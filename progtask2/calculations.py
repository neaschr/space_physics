import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors 
import datareader
from optical_depth import optical_depth
from euvflux import EUV_flux
from irradiance import irradiance_at_wavelength

if __name__ == '__main__':

    #Importing the data for absorption crossections
    names_of_columns = [
        'wavelength [m]', 
        'absorption cross section N2 [m^2]',
        'absorption cross section O  [m^2]', 
        'absorption cross section O2 [m^2]'
    ]

    #Finding the absorption crossesctions
    absorp_crosssec_data = datareader.data_reader('phot_abs.dat', names_of_columns)
    absorption_crosssecs = absorp_crosssec_data.iloc[:, 1:].to_numpy().T

    #Finding the wavelengths
    wavelengths = absorp_crosssec_data['wavelength [m]'].to_numpy()

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

    #Having all species number densities as one array
    N_all_species = number_densities = np.vstack([N_atom_oxy, N_molec_nitro, N_molec_oxy])
    
    #Defining array for the heights, given in unit meters
    height_in_meters = data['Height_km'].to_numpy() * 1e3
    z_0_values = height_in_meters[100:]

    #Making an array containing the different solar zenith angles we want to use
    sol_zen_ang = np.arange(0, 76, 15)
    sol_zen_ang = np.append(sol_zen_ang, 85)
    
    #Looping over the angles and plotting
    for angle in sol_zen_ang:

        #Finding optical depth for each z_0 and wavelength
        optical_depth_matrix = np.array([
        [optical_depth(z_0, angle, absorption_crosssecs[:, i], N_all_species, height_in_meters) for i in range(len(wavelengths))]
        for z_0 in z_0_values])

        #Use LogNorm for logarithmic color scaling
        norm = mcolors.LogNorm(vmin=0.01, vmax= 10**4) 

        #Plotting
        plt.figure(figsize=(10, 6))
        plt.pcolormesh(wavelengths * 1e9, z_0_values / 1e3, optical_depth_matrix, cmap='inferno', norm=norm)
        plt.colorbar(label='Optical Depth')
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Altitude (km)')
        plt.title(f'Optical Depth with $\\chi = {angle}^\\circ$')
        plt.show()

        #Defining column names and reading data
        colums = ['Time (yyyyDDD)', 'Wavelength (nm)', 'Irradiance (W/m^2/nm)', 'Uncertainty']
        irr_data = datareader.data_reader2(datafile_name = 'fism_daily_hr19990216.dat', column_names = colums)
        
        #Finding the data for wavelength and for irradiance. Converting to SI units
        wavelengths_m = irr_data['Wavelength (nm)'].to_numpy() * 1e9
        irradiances = irr_data['Irradiance (W/m^2/nm)'].to_numpy() * 1e9

        #Finding the flux, using both EUV_flux function and irradiance_at_wavelength function
        photon_flux_matrix = np.array([
            [EUV_flux(z_0, angle, absorption_crosssecs[:, i], N_all_species, height_in_meters, 
                      irradiance_at_wavelength(wavelengths[i], wavelengths_m, irradiances), 
                      wavelengths[i]) for i in range(len(wavelengths))] for z_0 in z_0_values])

        #Plotting
        plt.figure(figsize=(10, 6))
        plt.pcolormesh(wavelengths * 1e9, z_0_values / 1e3, photon_flux_matrix, cmap='inferno')
        plt.colorbar(label='Photon flux')
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Altitude (km)')
        plt.title(f'Solar EUV-flux with $\\chi = {angle}^\\circ$')
        plt.show()        
        