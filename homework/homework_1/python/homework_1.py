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
print(f'The fractional density change is {round(drho_rho_c*100,3)} %')

# Problem 1d
u_d = 980
drho_rho_d = -rho/(gamma*p)*u_d**2*du_u
print(f'The fractional density change is {round(drho_rho_d*100,3)} %')

magnitude = drho_rho_d/drho_rho_c
print(f'The fractional density change in part (d) is {round(magnitude,0)} times that of part (c)')

# Problem 3b
p_init = 1380 # kPa
T_init = 500 # K
gamma = 1.4
R_air = 287 # J/Kg K
# Assuming ideal gas (inherent b/c Isentropic), can use equation of state to get
# all initial conditions
rho_init = p_init*1000/(T_init*R_air)
print(f'Initial conditions:\n\tPressure = {p_init} kPa\n\tTemperature = {T_init} K\n\tDensity = {rho_init} kg/m^3')

# Given initial conditions and assuming an isentropic process,
# can assume that P/\rho^{\gammma} is constant.

isentropic_constant = p_init/rho_init**gamma
print(isentropic_constant)

tunnel_data = pd.read_csv('../HW_1_data.csv')
times = tunnel_data["Test Time (ms)"]
pressures = tunnel_data["Test-Section Pressure (Pa)"]
times = times.values.tolist()
pressures = pressures.values.tolist()


densities = [(pressure/isentropic_constant)**(1/gamma) for pressure in pressures]
temperatures = [pressure*1000/(density*R_air) for pressure,density in zip(pressures,densities)]

# TODO: Update time for plots

fig, ax = plt.subplots()
ax.plot(times,pressures,'.')
ax.set_xlabel('Time [ms]')
ax.set_ylabel('Pressure [kPa]')
plt.title('Pressure vs. Time')
plt.show()
plt.close()

fig, ax = plt.subplots()
ax.plot(times,densities,'.')
ax.set_xlabel('Time [ms]')
ax.set_ylabel(r'Density $[\frac{kg}{m^3}$]')
plt.title('Density vs. Time')
plt.show()
plt.close()

fig,ax = plt.subplots()
ax.plot(times,temperatures,'.')
ax.set_xlabel('Time [ms]')
ax.set_ylabel('Temperature [K]')
plt.title('Temperature vs. Time')
plt.show()
plt.close()

# Problem 3c
useful_pressures = [pressures[idx] for idx,time in enumerate(times) if (time >= 25) and (time <=75)]
mean = statistics.mean(useful_pressures)
deviation = statistics.stdev(useful_pressures)
print(f'The standard deviation of pressure about the mean between',
    f't = 25 ms and t = 75 ms is {round(deviation,3)} kPa')
print(f'The mean pressure between t = 25 ms and t = 75 ms is {round(mean,3)} kPa')

# Problem 3d
# At t= 50 ms, what is the entropy difference 
# between the upstream and test-section air? 
cp = 1000 # J/kg K
R = 287 # J/kg K
T_1 = T_init
p_1 = p_init
closest_time = min(times, key=lambda x:abs(x-50))
idx_50 = [idx for idx,time in enumerate(times) if time == closest_time]
print(idx_50)
T_2 = temperatures[idx_50[0]]
p_2 = pressures[idx_50[0]]
print(T_2,p_2)
ds_50 = cp*np.log(T_2/T_1) - R*np.log(p_2/p_1)
print(f'Change in entropy = {ds_50}')

print(f'T_1 = {T_1}, T_2 = {T_2}, p_2 = {p_1}, p_2 = {p_2} ')
print(densities[idx_50[0]]*R_air*T_2)
print(p_2*1000)

# Problem 3e
# For p_init = 1380 kPa, calculate the upstream starting temperature for
# which the air in the test section is equal to the liquefaction temperature
# at t = 50 ms
T_liq = 79 # http://hyperphysics.phy-astr.gsu.edu/hbase/thermo/liqair.html#:~:text=If%20air%20is%20cooled%20at,liquification%20point%20of%20pure%20nitrogen.
T_upstream = T_liq/(p_2/p_1)**((gamma-1)/gamma)
print(T_upstream)

# Problem 3f
Ma = 6
temp_40km = 250 # K