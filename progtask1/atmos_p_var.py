from datareading import data_reader
from altvar import altitude_var_g
from scipy import constants as spc
from matplotlib import pyplot as plt
import numpy as np

#We will calculate the atmospherc pressure variation 
def barometric_equation(pressure_0: float, mass: np.ndarray, heights: np.ndarray, grav_acc: np.ndarray, temperatures: np.ndarray) -> np.ndarray:
    '''Calculates atmospheric pressure, returns pressure array with units of pascal''' 

    k_B = spc.k
    delta_height = heights[1] - heights[0]
    delta_height = delta_height * 1e3 #converting to meters
    pressure = np.zeros(len(heights))
    pressure[0] = pressure_0 #Pa

    for n in range(1, len(heights)):
        exponent = - (mass[n] * grav_acc[n]) / (k_B * temperatures[n]) * delta_height
        pressure[n] = pressure[n-1] * np.exp(exponent)
        # pressure = pressure_0 * np.exp(- np.sum (((mass * grav_acc) / (k_B * temperatures)) * delta_height))

    return pressure

#Atomic mass unit
atomic_mass_unit = spc.physical_constants['atomic mass constant'][0] #kg

#reading the data file containing data about stuff
data = data_reader('MSIS.dat')

#Finding the number density as a sum of each species number density, given as number per cm cubed
number_density = (data['Oksygen-atoms_cm-3'].to_numpy() + data['Nitrogen-molecules_cm-3'].to_numpy() 
+ data['Oksygen-molecules_cm-3'].to_numpy())

#Calculating the average particle mass using the mass density and number density, mass in grams
avg_particle_mass = data['Mass_density_g_per_cm3'].to_numpy() / number_density

#Converting average particle mass to kilograms
avg_particle_mass = avg_particle_mass * 10**(-3)

#Temperatures
temperature_K = data['Temperature_neutral_K'].to_numpy()

#Using the function for gravitational acceleration from altvar.py. Converting the heights from km to m
grav_acc_var = altitude_var_g(heights= data['Height_km']*1000)

#Converting the heights to a numpy array
heights_in_km = data['Height_km'].to_numpy()

#Defining the pressure at z = 0
pressure_at_0 = 101300 #Pa

#Finding the pressure array for the atmosphere
atmospheric_pressure = barometric_equation(pressure_at_0, avg_particle_mass, heights_in_km, grav_acc_var, temperature_K)

#We will also calculate the pressure according to the ideal gas law using the temperatures and density profiles in the data file
def ideal_gas_law(number_density_per_cm3: np.ndarray, temperatures: np.ndarray) -> np.ndarray:

    k_B = spc.k
    number_density_per_meter3 = number_density_per_cm3 * 1e6
    pressure = (number_density_per_meter3 * k_B * temperatures) 

    return pressure

atmospheric_pressure_ideal = ideal_gas_law(number_density, temperature_K)

#Plotting...
fig, ax = plt.subplots(figsize = (9, 7))
ax.semilogx(atmospheric_pressure, heights_in_km, linestyle = 'solid', label = 'Barometric equation')
ax.semilogx(atmospheric_pressure_ideal, heights_in_km, linestyle = 'dashed', label = 'Ideal gas law')

#Plot aesthetics
ax.set_ylabel("Height (km)")
ax.set_xlabel("Pressure (Pa)")
ax.set_title("Pressure as a function of height")
ax.grid(True)

plt.legend()
plt.show()