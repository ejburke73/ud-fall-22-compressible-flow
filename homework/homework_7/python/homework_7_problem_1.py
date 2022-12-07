# Compressible Flow
# AEE 553
# Homework 7 - Problem 1
# Evan Burke

import numpy as np
from matplotlib import pyplot as plt
import shocks as ns
import oblique as os
import isentropic as isen
from scipy.optimize import fsolve
from scipy.optimize import root
gamma = 1.4
R = 287
d_tube = 9.75 # in
l_tube = 82 # ft
d_nozzle = 9.75 # in
l_nozzle = 117 # in
p0 = 4000000 # pa
T0 = 500 # K
Tb = 295 # K
pb = 2533.45 # Pa, from HW6
a1 = isen.get_sonic_velocity(T=Tb)
a4 = isen.get_sonic_velocity(T=T0)

# a
def shock_tube_pressures(p2_p1=None,a1=None,a4=None,gamma=1.4):
    p4_p1 = p2_p1 * (1 - (gamma-1)*(a1/a4)*(p2_p1-1) / (2*gamma*(2*gamma + (gamma+1)*(p2_p1-1)))**0.5)**((-2*gamma)/(gamma-1))
    return p4_p1

p2_p1 = 17.975422520 # from MATLAB

T2_T1 = p2_p1 * (((gamma+1)/(gamma-1)+p2_p1) / (1 + (gamma+1)/(gamma-1)*p2_p1))
print(f'T2/T1 = {T2_T1}')
rho2_rho1 = ( (1 + (gamma+1)/(gamma-1)*p2_p1) / ((gamma+1)/(gamma-1)+p2_p1))
print(f'rho2/rho1 = {rho2_rho1}')
W = a1 * ((gamma+1)/(2*gamma) *(p2_p1-1)+1)**0.5
print(f'W = {W}')
up = a1/gamma * (p2_p1 -1) * (((2*gamma)/(gamma+1))/(p2_p1 + (gamma-1)/(gamma+1)))**0.5
print(f'up = {up}')
Ms = W/a1
print(f'Ms = {Ms}')

p2 = p2_p1 * pb
print(f'p2 = {p2}')

T2 = T2_T1*Tb
print(f'T2  {T2}')

a2 = isen.get_sonic_velocity(T=T2)
M2 = up/a2
print(f'M2 = {M2}')

p2t = isen.get_total_pressure(M=M2,p=p2)
p2t_model = ns.get_total_pressure_normal_shock(M1=M2,p1_t=p2t)

l_nozzle = l_nozzle * 0.0254
print(f'Nozzle length = {l_nozzle} m')

#where is the probe
time_to_nozzle_NS = l_nozzle/W
print(f'Time for NS to hit exit of nozzle = {time_to_nozzle_NS}')

time_to_nozzle_CS = l_nozzle/up
print(f'Time for CS to hit exit of nozzle = {time_to_nozzle_CS}')


# Expansion Wave
p0 = 1000000
pb = 101000
p2_p1 = 3.2855668652704 # from MATLAB

p3_p4 = p2_p1 * p0/pb
print(f'p3/p4 = {p3_p4}')
rho3_rho4 = p3_p4**(1/gamma)
print(f'rho3/rho4 = {rho3_rho4}')
T3_T4 = p3_p4**((gamma-1)/gamma)
print(f'T3/T4 = {T3_T4}')
