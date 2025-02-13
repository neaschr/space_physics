from datareading import data_reader
from scaleheights import scale_height
from altvar import altitude_var_g
from scipy import constants as spc
from matplotlib import pyplot as plt
import numpy as np

#I will calculate the scale height for molecular O, molecular N and atomic O at 120 and 600 km altitude

#Approximated gravitational acceleration
g_approx = 9.28 #m/s^2

#Atomic mass unit
atomic_mass_unit = spc.physical_constants['atomic mass constant'][0] #kg

#Defining the mass of each atom or molecule. This is for individual scale heights
atomic_O_mass_kg = 16 * atomic_mass_unit
molecular_O_mass_kg = 32 * atomic_mass_unit
molecular_N_mass_kg = 28 * atomic_mass_unit

#reading the data file 
data = data_reader('MSIS.dat')

#Finding the number density as a sum of each species number density, given as number per cm cubed
number_density = data['Oksygen-atoms_cm-3'].to_numpy() + data['Nitrogen-molecules_cm-3'].to_numpy() 
+ data['Oksygen-molecules_cm-3'].to_numpy()

#Calculating the average particle mass using the mass density and number density, mass in grams
avg_particle_mass = data['Mass_density_g_per_cm3'].to_numpy() / number_density

#Converting average particle mass to kilograms
avg_particle_mass = avg_particle_mass * 10**(-3)

#Temperatures
temperature_K = data['Temperature_neutral_K'].to_numpy()

#Finding the scale height at 120 and 600 km of altitude using the approximated value of g = 9.28
scaleheigts_calculated_1 = scale_height(temperature= temperature_K, particle_mass= avg_particle_mass, grav_acc= g_approx)
H_at_120km_1 = scaleheigts_calculated_1[120]
H_at_600km_1 = scaleheigts_calculated_1[600]

#Using the function for gravitational acceleration from altvar.py. Converting the heights from km to m
grav_acc_var = altitude_var_g(heigths= data['Height_km']*1000)

#Using the altitude varying gravitational acceleration when calculating scale height
scaleheigts_calculated_2 = scale_height(temperature= temperature_K, particle_mass= avg_particle_mass, grav_acc= grav_acc_var)
H_at_120km_2 = scaleheigts_calculated_2[120]
H_at_600km_2 = scaleheigts_calculated_2[600]

print(f'Scale height at 120 km using g=9.28: H = {H_at_120km_1}')
print(f'Scale height at 600 km using g=9.28: H = {H_at_600km_1}')
print(f'Scale height at 120 km using varying g: H = {H_at_120km_2}')
print(f'Scale height at 600 km using varying g: H = {H_at_600km_2}')

#We see that as altitude increases, the difference in calculated scale height between the approximations for g also increase

#Converting the heights to a numpy array
heights_in_km = data['Height_km'].to_numpy()

#Plotting for fun
fig, ax = plt.subplots(figsize = (9, 7))
ax.plot(scaleheigts_calculated_1, heights_in_km, linestyle = 'dashed', label = r'$g = 9.28 ms^{-2}$')
ax.plot(scaleheigts_calculated_2, heights_in_km, linestyle = 'solid', label = r'$g(z) = g(0) \frac{R_E^{2}}{(R_E + z)^2}$')

#Plot aesthetics
ax.set_ylabel("Height (km)")
ax.set_xlabel("Scale height (km)")
ax.set_title("Scale height as a function of height")
ax.grid(True)

plt.legend()
plt.show()
