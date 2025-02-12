from matplotlib import pyplot as plt
from datareading import data_reader
from scaleheights import scale_height

#using the data reading function to read the data and return a dataframe with the correct collumn names
data = data_reader('MSIS.dat')

#Making arrays of the data i will use
temperatures = data['Temperature_neutral_K'].to_numpy() #Temperatures in Kelvin

#Finding the number density as a sum of each species number density, given as number per cm cubed
number_density = data['Oksygen-atoms_cm-3'].to_numpy() + data['Nitrogen-molecules_cm-3'].to_numpy() 
+ data['Oksygen-molecules_cm-3'].to_numpy()

#Calculating the average particle mass using the mass density and number density, mass in grams
avg_particle_mass = data['Mass_density_g_per_cm-3'].to_numpy() / number_density

#Converting average particle mass to kilograms
avg_particle_mass = avg_particle_mass * 10**(-3)

#Finding scale height for each 1 km interval
calculated_scaleheights = scale_height(temperature=temperatures, particle_mass=avg_particle_mass)

#Converting the heights to a numpy array
heights_in_km = data['Height_km'].to_numpy()

#Plotting
fig, ax = plt.subplots(figsize = (9, 6))
ax.plot(calculated_scaleheights, heights_in_km)

#Plot aesthetics
ax.set_ylabel("Height (km)")
ax.set_xlabel("Scale height (km)")
ax.set_title("Scale height as a function of height")
ax.grid(True)

# plt.legend()
plt.show()

#Finding scale height at ground level
scale_height_at_0 = calculated_scaleheights[0]

#Finding scale height at 100km of altitude
scale_height_at_100 = calculated_scaleheights[100]

#Printing the information
print(f'Scale height at 0 m above ground: H = {scale_height_at_0}km \n Scale height at 100km above ground: H = {scale_height_at_100} km')