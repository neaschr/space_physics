from datareading import data_reader
from scaleheights import scale_height
from scaleheights import scale_height_experimental
from altvar import altitude_var_g
from scipy import constants as spc
from matplotlib import pyplot as plt
import numpy as np

#We will compare different approximations for scale height for different species of particle

#Atomic mass unit
atomic_mass_unit = spc.physical_constants['atomic mass constant'][0] #kg

#Defining the mass of each atom or molecule. This is for individual scale heights
atomic_O_mass_kg = 16 * atomic_mass_unit
molecular_O_mass_kg = 32 * atomic_mass_unit
molecular_N_mass_kg = 28 * atomic_mass_unit

#reading the data file 
data = data_reader('MSIS.dat')

#Approximated gravitational acceleration
g_approx = 9.28 #m/s^2

#Using the function for gravitational acceleration from altvar.py. Converting the heights from km to m
grav_acc_var = altitude_var_g(heigths= data['Height_km']*1000)

#Converting the heights to a numpy array
heights_in_km = data['Height_km'].to_numpy()

#importing the temperature
temperatures_K = data['Temperature_neutral_K'].to_numpy()

#Importing the number densities for each particle
number_dens_molecular_N = data['Nitrogen-molecules_cm-3'].to_numpy()
number_dens_atomic_O = data['Oksygen-atoms_cm-3'].to_numpy()
number_dens_molecular_O = data['Oksygen-molecules_cm-3'].to_numpy()

#Calculating the scale height for each individual particle species using the appproximated g
scale_H_N2_g_approx = scale_height(temperatures_K, molecular_N_mass_kg, grav_acc = g_approx)
scale_H_O_g_approx = scale_height(temperatures_K, atomic_O_mass_kg, grav_acc = g_approx)
scale_H_O2_g_approx = scale_height(temperatures_K, molecular_O_mass_kg, grav_acc = g_approx)

#Calculating the scale height for each individual particle species using the altitude varying g
scale_H_N2_altvar = scale_height(temperatures_K, molecular_N_mass_kg, grav_acc = grav_acc_var)
scale_H_O_altvar = scale_height(temperatures_K, atomic_O_mass_kg, grav_acc = grav_acc_var)
scale_H_O2_altvar = scale_height(temperatures_K, molecular_O_mass_kg, grav_acc = grav_acc_var)

#Calculating scale height for each individual particle species using experimental scale height
scaleH_N2_experimental = scale_height_experimental(number_density= number_dens_molecular_N)
scaleH_O2_experimental = scale_height_experimental(number_density= number_dens_molecular_O)
scaleH_O_experimental = scale_height_experimental(number_density= number_dens_atomic_O)

#Plotting for molecular nitrogen
fig, ax = plt.subplots(figsize = (9, 7))
ax.plot(scale_H_N2_g_approx, heights_in_km, linestyle = 'dashed', label = r'$g = 9.28 ms^{-2}$')
ax.plot(scale_H_N2_altvar, heights_in_km, linestyle = 'solid', label = r'$g(z) = g(0) \frac{R_E^{2}}{(R_E + z)^2}$')
ax.plot(scaleH_N2_experimental, heights_in_km, linestyle = 'dotted', label = r'$H = - \frac{n_{N2}}{\partial n_{N2} / \partial z}$')

#Plot aesthetics
ax.set_ylabel("Height (km)")
ax.set_xlabel("Scale height for N2 (km)")
ax.set_title("Scale height for molecular nitrogen as a function of height")
ax.grid(True)

plt.legend()
plt.show()

#Plotting for molecular oxygen
fig, ax = plt.subplots(figsize = (9, 7))
ax.plot(scale_H_O2_g_approx, heights_in_km, linestyle = 'dashed', label = r'$g = 9.28 ms^{-2}$')
ax.plot(scale_H_O2_altvar, heights_in_km, linestyle = 'solid', label = r'$g(z) = g(0) \frac{R_E^{2}}{(R_E + z)^2}$')
ax.plot(scaleH_O2_experimental, heights_in_km, linestyle = 'dotted', label = r'$H = - \frac{n_{O2}}{\partial n_{O2} / \partial z}$')

#Plot aesthetics
ax.set_ylabel("Height (km)")
ax.set_xlabel("Scale height for O2 (km)")
ax.set_title("Scale height for molecular oxygen as a function of height")
ax.grid(True)

plt.legend()
plt.show()

#Plotting for atomic oxygen
fig, ax = plt.subplots(figsize = (9, 7))
ax.plot(scale_H_O_g_approx, heights_in_km, linestyle = 'dashed', label = r'$g = 9.28 ms^{-2}$')
ax.plot(scale_H_O_altvar, heights_in_km, linestyle = 'solid', label = r'$g(z) = g(0) \frac{R_E^{2}}{(R_E + z)^2}$')
ax.plot(scaleH_O_experimental, heights_in_km, linestyle = 'dotted', label = r'$H = - \frac{n_O}{\partial n_O / \partial z}$')

#Plot aesthetics
ax.set_ylabel("Height (km)")
ax.set_xlabel("Scale height for O (km)")
ax.set_title("Scale height for atomic oxygen as a function of height")
ax.grid(True)

plt.legend()
plt.show()