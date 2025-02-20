from datareading import data_reader
from scaleheights import scale_height
from altvar import altitude_var_g
from scipy import constants as spc
from matplotlib import pyplot as plt

#Atomic mass unit
atomic_mass_unit = spc.physical_constants['atomic mass constant'][0] #kg

#reading the data file containing data about stuff
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

#Using the function for gravitational acceleration from altvar.py. Converting the heights from km to m
grav_acc_var = altitude_var_g(heigths= data['Height_km']*1000)

#Using the altitude varying gravitational acceleration when calculating scale height
scaleheigts_calculated = scale_height(temperature= temperature_K, particle_mass= avg_particle_mass, grav_acc= grav_acc_var)

#Converting the heights to a numpy array
heights_in_km = data['Height_km'].to_numpy()

#Plotting for fun
fig, ax = plt.subplots(figsize = (9, 7))
ax.plot(scaleheigts_calculated, heights_in_km, linestyle = 'solid')


#Plot aesthetics
ax.set_ylabel("Height (km)")
ax.set_xlabel("Scale height (km)")
ax.set_title("Scale height as a function of height")
ax.grid(True)

plt.show()


