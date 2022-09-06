#!/usr/bin/env python3
from operator import index
import pandas as pd
from matplotlib import pyplot as plt
import statistics
import numpy as np

rho = 1.23
p = 1.01e5
du_u = 0.01
gamma = 1.4

# Problem 1c
u_c = 63
drho_rho_c = -rho/(gamma*p)*u_c**2*du_u
print('Problem 1c:')
print(f'The fractional density change is {round(drho_rho_c*100,3)} %\n')

# Problem 1d
u_d = 980
drho_rho_d = -rho/(gamma*p)*u_d**2*du_u
print('Problem 1d:')
print(f'The fractional density change is {round(drho_rho_d*100,3)} %')

magnitude = drho_rho_d/drho_rho_c
print(f'The fractional density change in part (d) is {round(magnitude,0)} times that of part (c)\n')

# Problem 3b
p_init = 1380000 # Pa
T_init = 500 # K
gamma = 1.4
R_air = 287 # J/Kg K

# Assuming ideal gas (inherent b/c Isentropic), can use equation of state to get
# all initial conditions
rho_init = p_init/(T_init*R_air) # Convert kPa to Pa for dimension analysis

print('Problem 3b:')
print(f'Initial conditions:\n\tPressure = {p_init} Pa\n\tTemperature = {T_init} K\n\tDensity = {round(rho_init,2)} kg/m^3')

# Given initial conditions and assuming an isentropic process,
# can assume that P/\rho^{\gammma} is constant.

isentropic_constant = p_init/rho_init**gamma
print(f'Isentropic constant = {isentropic_constant}\n')

tunnel_data = pd.read_csv('../HW_1_data.csv')
times = tunnel_data["Test Time (ms)"]
pressures = tunnel_data["Test-Section Pressure (Pa)"] 
times = times.values.tolist()
pressures = pressures.values.tolist()

densities = [(pressure/isentropic_constant)**(1/gamma) for pressure in pressures]
temperatures = [pressure/(density*R_air) for pressure,density in zip(pressures,densities)]

fig, ax = plt.subplots()
ax.plot(times,pressures,'.')
ax.set_xlabel('Time [ms]')
ax.set_ylabel('Pressure [Pa]')
plt.title('AFRL Ludwieg Tube\nPressure vs. Time')
plt.savefig('../images/pressure_vs_time.png')
plt.close()

fig, ax = plt.subplots()
ax.plot(times,densities,'.')
ax.set_xlabel('Time [ms]')
ax.set_ylabel(r'Density $[\frac{kg}{m^3}$]')
plt.title('AFRL Ludwieg Tube\nDensity vs. Time')
plt.savefig('../images/density_vs_time.png')
plt.close()

fig,ax = plt.subplots()
ax.plot(times,temperatures,'.')
ax.set_xlabel('Time [ms]')
ax.set_ylabel('Temperature [K]')
plt.title('AFRL Ludwieg Tube\nTemperature vs. Time')
plt.savefig('../images/temperature_vs_time.png')
plt.close()


# Problem 3c
useful_pressures = [pressures[idx] for idx,time in enumerate(times) if (time >= 25) and (time <=75)]
mean = statistics.mean(useful_pressures)
deviation = statistics.stdev(useful_pressures)
print('Problem 3c:')
print(f'The standard deviation of pressure about the mean between',
    f't = 25 ms and t = 75 ms is {round(deviation,3)} Pa')
print(f'The mean pressure between t = 25 ms and t = 75 ms is {round(mean,2)} Pa\n')

# Problem 3d
# At t= 50 ms, what is the entropy difference 
# between the upstream and test-section air? 
cp = 1000 # J/kg K
R = 287 # J/kg K
T_1 = T_init
p_1 = p_init
closest_time = min(times, key=lambda x:abs(x-50))
idx_50 = [idx for idx,time in enumerate(times) if time == closest_time]
T_2 = temperatures[idx_50[0]]
p_2 = pressures[idx_50[0]]
rho_2 = densities[idx_50[0]]
ds_50 = cp*np.log(T_2/T_1) - R*np.log(p_2/p_1)
print('Problem 3d:')
print(f'T_1 = {round(T_1,4)}, T_2 = {round(T_2,2)}, p_1 = {round(p_1, 2)}, p_2 = {round(p_2,2)} ')
print(f'Change in entropy = {round(ds_50,2)} J/kg K\n')


# Problem 3e
# For p_init = 1380 kPa, calculate the upstream starting temperature for
# which the air in the test section is equal to the liquefaction temperature
# at t = 50 ms
# Assuming that the pressure at t = 50 ms is still the temperature from the .xls file
print('Problem 3e:')
T_liq = 77 # http://www.edubilla.com/invention/liquefaction-of-air/
T_upstream = T_liq/(p_2/p_1)**((gamma-1)/gamma)
print(f'The upstream temperature required to cause liquefaction temperature at t = 50 ms is {round(T_upstream,2)} K')

